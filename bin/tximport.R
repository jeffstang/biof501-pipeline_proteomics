#!/usr/bin/env Rscript

library(tximport)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
quant_dirs <- unlist(strsplit(args[1], split=" "))
tx2gene_input <- args[2] 
output_prefix <- args[3]

# Check if input files exist
if (!file.exists(metadata_file)) stop("Metadata file does not exist")
if (!file.exists(tx2gene_input)) stop("tx2gene file does not exist")

print("Reading metadata and tx2gene files...")

# Read metadata
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Read tx2gene file
# This file will handle the naming conversion of genCode transcript IDs to gene symbols
tx2gene <- read.delim(tx2gene_input, header = FALSE, stringsAsFactors = FALSE, sep = "")
colnames(tx2gene) <- c("TXNAME", "GENEID")

# Creating a print statement to show that 
print("Processing sample information...")

# This will list out the expected sample outDirs generated from SALMON_QUANT
# For this demo, we will expect 4 runs so 4 directories will be listed.
# I want to extract the gene-level counts files of each sample
quant_genes_files <- unlist(lapply(quant_dirs, function(d) {
  f <- list.files(d, pattern = "quant.genes.sf", full.names = TRUE)
  names(f) <- d
  return(f)
}))

print("Importing quantifications with tximport...")

# Import quantifications
# I followed the tximport tutorial vignette on creating a txi object
# Because limma-voom does not use the offset matrix stored in y$offset, 
# They recommend using the scaled counts generated from abundances
# I used lengthscaledTPM here:
txi <- tximport(quant_genes_files, type = "salmon", tx2gene = tx2gene, txOut = TRUE, dropInfReps=TRUE,countsFromAbundance = "lengthScaledTPM")

# 1) Expected output will be a gene x sample counts matrix
#    Rownames : Gene symbols
#    Colnames : Sample by SRR IDs 
# 2) Save the TXI object into an .RDS
# Although I may only use txi$counts for the downstream analysis, I chose to provide
# flexibility for future runs 
saveRDS(txi, file = paste0(output_prefix, "_txi.rds"))

# print("Process completed successfully.")