#!/usr/bin/env Rscript

library(tximport)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
quant_dirs <- args[1]
metadata_file <- args[2]
tx2gene_input <- args[3] 
output_prefix <- args[4]

# Check if input files exist
if (!file.exists(metadata_file)) stop("Metadata file does not exist")
if (!file.exists(tx2gene_input)) stop("tx2gene file does not exist")

print("Reading metadata and tx2gene files...")
# Read metadata
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
# Read tx2gene file
tx2gene <- read.delim(tx2gene_input, header = FALSE, stringsAsFactors = FALSE)
colnames(tx2gene) <- c("TXNAME", "GENEID")

print("Processing sample information...")
# Read sample information
samples <- list.dirs(quant_dirs, full.names = TRUE, recursive = FALSE)
files <- file.path(samples, "quant.genes.sf")
names(files) <- basename(samples)

# Check if all quant.genes.sf files exist
if (!all(file.exists(files))) stop("Some quant.genes.sf files are missing")

# Match sample names with metadata
sample_names <- basename(samples)
metadata_subset <- metadata[metadata$SRR_ID %in% sample_names, ]
metadata_subset <- metadata_subset[match(sample_names, metadata_subset$SRR_ID), ]

if (nrow(metadata_subset) != length(samples)) {
  stop("Mismatch between number of samples and metadata entries")
}

# Rename files to match metadata
names(files) <- metadata_subset$geo_accession

print("Importing quantifications with tximport...")
# Import quantifications
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = FALSE)

# Add metadata to the txi object
txi$metadata <- metadata_subset

print("Saving txi object...")
# Save txi object for later use
saveRDS(txi, file = paste0(output_prefix, "_txi.rds"))

print("Process completed successfully.")