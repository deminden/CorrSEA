suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(parallel)
  library(recount3)
  library(SummarizedExperiment)
})

# Allow long TCGA downloads to complete.
options(timeout = 3600)

# Set folder paths
data_folder <- "data"
output_folder_norm <- file.path(data_folder, "TCGA_normalised_counts_cancers")
output_folder_vst <- file.path(data_folder, "TCGA_vst_counts_cancers")
annotation <- "gencode_v29"
recount3_url <- getOption(
  "recount3_url",
  "https://recount-opendata.s3.amazonaws.com/recount3/release"
)

parse_args <- function(default_workers = 1) {
  args <- commandArgs(trailingOnly = TRUE)
  workers_arg <- NULL
  cancers_arg <- NULL

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (grepl("^--workers=", arg)) {
      workers_arg <- sub("^--workers=", "", arg)
      i <- i + 1
      next
    }

    if (identical(arg, "--workers")) {
      if (i == length(args)) {
        stop("--workers must be followed by a positive integer.")
      }
      workers_arg <- args[[i + 1]]
      i <- i + 2
      next
    }

    if (grepl("^--cancers=", arg)) {
      cancers_arg <- sub("^--cancers=", "", arg)
      i <- i + 1
      next
    }

    if (identical(arg, "--cancers")) {
      if (i == length(args)) {
        stop("--cancers must be followed by one or more comma-separated cancer IDs.")
      }
      cancers_arg <- args[[i + 1]]
      i <- i + 2
      next
    }

    if (!startsWith(arg, "--") && is.null(workers_arg)) {
      workers_arg <- arg
      i <- i + 1
      next
    }

    stop(sprintf("Unknown argument: %s", arg))
  }

  if (is.null(workers_arg)) {
    workers <- default_workers
  } else {
    workers <- suppressWarnings(as.integer(workers_arg))
    if (is.na(workers) || workers < 1) {
      stop("Worker count must be a positive integer.")
    }
  }

  cancers <- NULL
  if (!is.null(cancers_arg)) {
    cancers <- trimws(unlist(strsplit(cancers_arg, ",")))
    cancers <- cancers[nzchar(cancers)]
    if (length(cancers) == 0) {
      stop("--cancers must contain at least one non-empty cancer ID.")
    }
  }

  list(workers = workers, cancers = cancers)
}

args <- parse_args()
workers <- args$workers
requested_cancers <- args$cancers

# Create output directories if they don't exist
dir.create(data_folder, showWarnings = FALSE)
dir.create(output_folder_norm, showWarnings = FALSE)
dir.create(output_folder_vst, showWarnings = FALSE)

# Confirm requested annotation is available.
available_annotations <- annotation_options("human")
if (!(annotation %in% available_annotations)) {
  stop(sprintf(
    "Annotation '%s' is not available. Available human annotations: %s",
    annotation,
    paste(available_annotations, collapse = ", ")
  ))
}

# Load input data
human_projects <- available_projects(recount3_url = recount3_url)
tcga_projects <- subset(
  human_projects,
  organism == "human" &
    file_source == "tcga" &
    project_type == "data_sources"
)

if (nrow(tcga_projects) == 0) {
  stop("No TCGA data source projects found in recount3.")
}

if (!is.null(requested_cancers)) {
  requested_cancers <- toupper(requested_cancers)
  missing_cancers <- setdiff(requested_cancers, tcga_projects$project)
  if (length(missing_cancers) > 0) {
    stop(sprintf(
      "Requested cancer IDs not found in recount3 TCGA projects: %s",
      paste(missing_cancers, collapse = ", ")
    ))
  }
  tcga_projects <- tcga_projects[tcga_projects$project %in% requested_cancers, ]
}

tcga_projects <- tcga_projects[order(tcga_projects$project), ]

sanitize_project_name <- function(project) {
  gsub("[^A-Za-z0-9_]+", "_", project)
}

# Function to normalize counts and save results
normalize_counts <- function(project_info) {
  cancer <- project_info$project
  cancer_file <- sanitize_project_name(cancer)

  # Define file paths
  norm_file <- file.path(
    output_folder_norm,
    paste0(cancer_file, "_", annotation, "_normalised_counts.tsv.gz")
  )
  vst_file <- file.path(
    output_folder_vst,
    paste0(cancer_file, "_", annotation, "_vst_counts.tsv.gz")
  )
  norm_needed <- !file.exists(norm_file)
  vst_needed <- !file.exists(vst_file)

  # Skip processing if requested output files already exist
  if (!norm_needed && !vst_needed) {
    cat(sprintf("Files for cancer '%s' already exist. Skipping processing.\n", cancer))
    return(list(status = "skipped", cancer = cancer, error = NA_character_))
  }

  result <- tryCatch({
    # Download recount3 data and convert coverage counts to read-style counts.
    rse <- create_rse(project_info, annotation = annotation, recount3_url = recount3_url)
    count_matrix <- round(transform_counts(rse))
    rownames(count_matrix) <- rownames(rse)

    # Extract colData
    sample_ids <- colnames(count_matrix)
    colData <- data.frame(row.names = sample_ids, condition = rep("tumour", length(sample_ids)))

    # Create DESeq2 dataset and estimate size factors for normalized counts.
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~1)
    dds <- estimateSizeFactors(dds)

    # Write normalized counts
    if (norm_needed) {
      norm_counts <- counts(dds, normalized = TRUE)
      norm_dt <- as.data.table(norm_counts, keep.rownames = "Ensembl_gene_ID")
      fwrite(norm_dt, file = norm_file, sep = "\t", compress = "gzip", row.names = FALSE)
    }

    # Perform VST transformation
    if (vst_needed) {
      vsd <- vst(dds)

      # Write VST-transformed data
      vst_counts <- assay(vsd)
      vst_dt <- as.data.table(vst_counts, keep.rownames = "Ensembl_gene_ID")
      fwrite(vst_dt, file = vst_file, sep = "\t", compress = "gzip", row.names = FALSE)
    }

    cat(sprintf("Processing complete for cancer '%s'.\n", cancer))
    list(status = "completed", cancer = cancer, error = NA_character_)
  }, error = function(e) {
    cat(sprintf("Processing failed for cancer '%s': %s\n", cancer, conditionMessage(e)))
    list(status = "failed", cancer = cancer, error = conditionMessage(e))
  })

  result
}

# Run normalization for each cancer in parallel
cat(sprintf("Using %d worker%s for TCGA cancer processing.\n",
            workers,
            ifelse(workers == 1, "", "s")))
cat(sprintf("Using recount3 annotation '%s'.\n", annotation))
cat(sprintf("Using recount3 URL '%s'.\n", recount3_url))

results <- mclapply(seq_len(nrow(tcga_projects)), function(i) {
  normalize_counts(tcga_projects[i, , drop = FALSE])
}, mc.cores = workers)

results_dt <- rbindlist(results)

completed_cancers <- results_dt[status == "completed", cancer]
skipped_cancers <- results_dt[status == "skipped", cancer]
failed_results <- results_dt[status == "failed"]

cat(sprintf("Completed cancers: %d\n", length(completed_cancers)))
cat(sprintf("Skipped cancers: %d\n", length(skipped_cancers)))
cat(sprintf("Failed cancers: %d\n", nrow(failed_results)))

if (nrow(failed_results) > 0) {
  cat("Failed cancers:\n")
  for (i in seq_len(nrow(failed_results))) {
    cat(sprintf("%s: %s\n", failed_results$cancer[[i]], failed_results$error[[i]]))
  }
}
