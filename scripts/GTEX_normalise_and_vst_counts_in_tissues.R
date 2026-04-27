suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(parallel)
})

# Allow long GTEx downloads to complete.
options(timeout = 3600)

# Set folder paths
data_folder <- "data"
output_folder_norm <- "Normalised_counts_tissues"
output_folder_vst <- "VST_counts_tissues"
vst_tissues <- character(0)

parse_worker_count <- function(default = 2) {
  args <- commandArgs(trailingOnly = TRUE)
  workers_arg <- NULL

  for (i in seq_along(args)) {
    if (grepl("^--workers=", args[[i]])) {
      workers_arg <- sub("^--workers=", "", args[[i]])
      break
    }

    if (identical(args[[i]], "--workers")) {
      if (i == length(args)) {
        stop("--workers must be followed by a positive integer.")
      }
      workers_arg <- args[[i + 1]]
      break
    }

    if (!startsWith(args[[i]], "--") && is.null(workers_arg)) {
      workers_arg <- args[[i]]
    }
  }

  if (is.null(workers_arg)) {
    return(default)
  }

  workers <- suppressWarnings(as.integer(workers_arg))
  if (is.na(workers) || workers < 1) {
    stop("Worker count must be a positive integer.")
  }

  workers
}

workers <- parse_worker_count()

# Create output directories if they don't exist
dir.create(data_folder, showWarnings = FALSE)
dir.create(output_folder_norm, showWarnings = FALSE)
if (length(vst_tissues) > 0) {
  dir.create(output_folder_vst, showWarnings = FALSE)
}

# Load input data
gene_reads_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.parquet"
annotation_url <- "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
gene_reads_file <- file.path(data_folder, basename(gene_reads_url))
annotation_file <- file.path(data_folder, basename(annotation_url))

download_if_missing <- function(url, file) {
  if (!file.exists(file)) {
    download.file(url, destfile = file, mode = "wb")
  }
}

download_if_missing(gene_reads_url, gene_reads_file)
download_if_missing(annotation_url, annotation_file)

# Read Parquet file
gene_reads_data <- arrow::read_parquet(gene_reads_file)
gene_reads_data <- as.data.frame(gene_reads_data, check.names = FALSE)

# Read annotation file
annotations <- fread(annotation_file)
annotations <- as.data.frame(annotations)

# Identify shared sample columns explicitly instead of relying on column position.
metadata_cols <- intersect(c("Name", "Description"), colnames(gene_reads_data))
sample_ids <- intersect(colnames(gene_reads_data), annotations$SAMPID)

if (length(sample_ids) == 0) {
  stop("No overlapping sample IDs found between the count matrix and the annotation file.")
}

annotations <- annotations[annotations$SAMPID %in% sample_ids, ]

# Split data by tissue
tissue_data <- split(annotations, annotations$SMTSD)

# Function to normalize counts and save results
normalize_counts <- function(tissue, tissue_df, gene_reads_data) {
  # Define file paths
  norm_file <- file.path(output_folder_norm, paste0(tissue, "_v11_normalised_counts.tsv.gz"))
  vst_file <- file.path(output_folder_vst, paste0(tissue, "_v11_vst_counts.tsv.gz"))
  make_vst <- tissue %in% vst_tissues
  norm_needed <- !file.exists(norm_file)
  vst_needed <- make_vst && !file.exists(vst_file)
  
  # Skip processing if requested output files already exist
  if (!norm_needed && !vst_needed) {
    cat(sprintf("Files for tissue '%s' already exist. Skipping processing.\n", tissue))
    return(NULL)
  }
  
  # Extract tissue-specific samples and counts
  tissue_samples <- tissue_df$SAMPID
  tissue_counts <- gene_reads_data[, c(metadata_cols, tissue_samples), drop = FALSE]
  
  # Extract count matrix and colData
  count_matrix <- as.matrix(tissue_counts[, tissue_samples, drop = FALSE])
  rownames(count_matrix) <- tissue_counts$Name
  colData <- data.frame(row.names = tissue_samples, condition = rep("control", length(tissue_samples)))
  
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
  
  cat(sprintf("Processing complete for tissue '%s'.\n", tissue))
}

# Run normalization for each tissue in parallel
cat(sprintf("Using %d worker%s for tissue processing.\n",
            workers,
            ifelse(workers == 1, "", "s")))

mclapply(names(tissue_data), function(tissue) {
  tissue_df <- tissue_data[[tissue]]
  normalize_counts(tissue, tissue_df, gene_reads_data)
}, mc.cores = workers)
