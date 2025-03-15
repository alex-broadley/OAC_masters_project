library(Seurat)
#devtools::install_github("miccec/yaGST")
#devtools::install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
library(dplyr)
library(scater)

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

#write.csv(tumoural_cells,"tumor_expression_mat.csv")

# Load CSV of epithelial counts matrix - genes = rows, columns = cells
normal_ids = read.csv("normal_cell_ids.csv")
normal_ids = normal_ids$index
normal_ids = gsub('-','.', normal_ids)

#epi_counts = read.csv("epi_raw_counts_matrix.csv")

# make gene-symbols row index - not a column
#rownames(epi_counts) = epi_counts$X
#epi_counts <- subset(epi_counts, select = -X)

epi_normal_counts = read.csv("epi_normal_raw_counts.csv")



# make gene-symbols row index - not a column
rownames(epi_normal_counts) = epi_normal_counts$X
epi_normal_counts <- subset(epi_normal_counts, select = -X)

epi_normal_counts = mutate_all(epi_normal_counts, function(x) as.numeric(as.character(x)))

#run full CNA pipeline
results_Ref <- pipelineCNA(epi_normal_counts,
                       norm_cell = normal_ids,
                       sample = 'tumor_class_w_normal',
                       SUBCLONES = FALSE,
                       ClonalCN = FALSE,
                       plotTree = FALSE
                       )

table(results_Ref$class)

#output to integrate in python
write.csv(results_Ref,"SCEVAN_results_wREF.csv", row.names = TRUE)
