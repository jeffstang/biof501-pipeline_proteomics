#!/usr/bin/env Rscript

library(tximport)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
quant_dirs <- unlist(strsplit(args[1], split=" "))
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
tx2gene <- read.delim(tx2gene_input, header = FALSE, stringsAsFactors = FALSE, sep = "")
colnames(tx2gene) <- c("TXNAME", "GENEID")

print("Processing sample information...")
# Read sample information
quant_genes_files <- unlist(lapply(quant_dirs, function(d) {
  f <- list.files(d, pattern = "quant.genes.sf", full.names = TRUE)
  names(f) <- d
  return(f)
}))

sample_names <- quant_dirs
names(quant_genes_files)
metadata_subset <- metadata[metadata$SRR_ID %in% sample_names, ]

print("Importing quantifications with tximport...")
# Import quantifications
txi <- tximport(quant_genes_files, type = "salmon", tx2gene = tx2gene, txOut = TRUE, dropInfReps=TRUE,countsFromAbundance = "lengthScaledTPM")

df <- data.frame(txi$abundance)

saveRDS(txi, file = paste0(output_prefix, "_txi.rds"))

# print("Process completed successfully.")