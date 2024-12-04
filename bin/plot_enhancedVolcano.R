#!/usr/bin/env Rscript

# Load Enhanced Volcano Package
# Load ggplot
library(EnhancedVolcano)
library(ggplot2)

# Command line input
args <- commandArgs(trailingOnly = TRUE)
top_table_csv <- args[1]

# Read in the top table given the input on the command line
top <- read.csv(top_table_csv, row.names=1)

# Create the EnhancedVolcano plot
volcano_plot <- EnhancedVolcano(top,
                lab = rownames(top),
                x = 'logFC',
                y = 'P.Value')

ggsave(filename="volcano_plot.png", plot=volcano_plot)