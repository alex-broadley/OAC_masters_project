library(Seurat)
#devtools::install_github("miccec/yaGST")
#devtools::install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
library(dplyr)
library(scater)

setwd("/Users/alex/Documents/University/Year 4/BIOL0041/SC_analysis")
write.csv(tumoural_cells,"tumor_expression_mata.csv")


# Load CSV of epithelial counts matrix - genes = rows, columns = cells
epi_counts = read.csv("epi_raw_counts_matrix.csv")

# make gene-symbols row index - not a column
rownames(epi_counts) = epi_counts$gene_symbols
epi_counts <- subset(epi_counts, select = -gene_symbols)

#run full CNA pipeline
results_noRef <- pipelineCNA(epi_counts,
                       SUBCLONES = TRUE,
                       ClonalCN = TRUE,
                       plotTree = TRUE
                       )

table(results_noRef$class)

#output to integrate in python
write.csv(results_noRef,"SCEVAN_results.csv", row.names = TRUE)