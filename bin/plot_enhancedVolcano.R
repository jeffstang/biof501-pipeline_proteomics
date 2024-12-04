#!/usr/bin/env Rscript

# Load Enhanced Volcano Package
library(EnhancedVolcano)

# Command line input
top_table_csv <- args[1]

# Read in the top table given the input on the command line
top <- read.csv(top_table_csv, row.names=1)

# Specify output file and resolution
png("volcano_plot.png", width = 800, height = 800, res = 300)

# Create the EnhancedVolcano plot
EnhancedVolcano(top,
                lab = rownames(top),
                x = 'logFC',
                y = 'P.Value')

# Close the graphical device
dev.off()