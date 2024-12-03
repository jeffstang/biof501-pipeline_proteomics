#!/usr/bin/env Rscript

library(tximport)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
quant_dirs <- unlist(strsplit(args[1], split=" "))
tx2gene_input <- args[2] 
output_prefix <- args[3]

# Check if input file exist
if (!file.exists(tx2gene_input)) stop("tx2gene file does not exist")

print("Reading tx2gene files...")

# This file will handle the naming conversion of genCode transcript IDs to gene symbols
tx2gene <- read.delim(tx2gene_input, header = FALSE, stringsAsFactors = FALSE, sep = "")
colnames(tx2gene) <- c("TXNAME", "GENEID")

print("Processing sample information...")

# This will list out the expected sample outDirs generated from SALMON_QUANT
# For this demo, we will expect 4 runs so 4 directories will be listed.
# Extract gene-level counts files of each sample
quant_genes_files <- unlist(lapply(quant_dirs, function(d) {
  f <- list.files(d, pattern = "quant.genes.sf", full.names = TRUE)
  names(f) <- d
  return(f)
}))

print("Importing quantifications with tximport...")

# Import quantifications
# I followed the tximport tutorial vignette on creating a txi object 
# They recommend using the scaled counts generated from abundances, however I do not do this
# As I am directly using the gene-level counts, which should already be scaled accordingly.
txi <- tximport(quant_genes_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE, countsFromAbundance = "no")

# 1) Expected output will be a gene x sample counts matrix
#    Rownames : Gene symbols
#    Colnames : Sample by SRR IDs 
# 2) Save the TXI object into an .RDS
# Although I may only use txi$counts for the downstream analysis, I chose to provide
# flexibility for future runs 
write.csv(txi$counts, file = paste0(output_prefix, "_txi.csv"))

# print("Process completed successfully.")_