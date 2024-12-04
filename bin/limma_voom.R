#!/usr/bin/env Rscript

# Load edgeR, this will also automatically load limma
library(edgeR)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
txi_file <- args[1]
metadata_file <- args[2]
out_prefix <- args[3]

# Check if input files exist:
if (!file.exists(txi_file)) stop("Salmon counts matrix does not exist")
if (!file.exists(metadata_file)) stop("Metadata file does not exist")

# Read in tximport object
txi <- read.csv(txi_file, row.names = 1, header = TRUE)

# Read in metadata of samples provided
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Make sure that the metadata and counts matrix follow the exact order 
# I made sure to order the counts matrix by the metadata
# For readability, it's just easier to keep track of conditions if they are grouped together
txi <- txi[, match(metadata$SRR_ID, colnames(txi))]

# Verify alignment
if (!all(metadata$SRR_ID == colnames(txi))) {
  stop("Sample names in metadata and counts matrix do not match after reordering!")
}

# Get the condition (usually some type of perturbation/treatment vs control)
condition <- factor(metadata$treatment)

# Store the counts into a DGEList so that it can be filtered on design
y <- DGEList(txi)

# Filter using the design information:
design <- model.matrix(~condition, data = metadata)
keep <- filterByExpr(y, design)
y <- y[keep, ]

# Normalization:
# Scale on TMM which has been found to perform well in comparative studies
# as stated in the limma User's Guide:
y <- calcNormFactors(y)

# Voom transformation:
# For troubleshooting, it is best practice to set plot=TRUE to look at the data
# but for the purpose of this pipeline, I did not declare this argument
v <- voom(y, design)

# Run limma:
fit <- lmFit(v, design)
fit <- eBayes(fit)
top <- topTable(fit, coef=ncol(design), number = nrow(fit))

write.csv(top, file = paste0(out_prefix, "_topTable.csv"))