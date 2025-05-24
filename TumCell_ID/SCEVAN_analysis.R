library(Seurat)
#devtools::install_github("miccec/yaGST")
#devtools::install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
library(dplyr)
library(scater)

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

# Load CSV of epithelial counts matrix - genes = rows, columns = cells
normal_ids = read.csv("normal_cell_ids.csv")
normal_ids = normal_ids$index
normal_ids = gsub('-','.', normal_ids)

#cehck have expected number of normal cells
normal_ids[which(normal_ids %in% colnames(epi_counts))]
length(normal_ids)
length(normal_ids[which(normal_ids %in% colnames(epi_normal_raw_counts))])
length(normal_ids[which(normal_ids %in% colnames(epi_und_normal_raw_counts))])

epi_und_counts = read.csv("epi_und_raw_counts.csv", row.names='X')
epi_und_TME_raw_counts = read.csv("epi_und_TME_raw_counts.csv", row.names = 'X')

epi_normal_raw_counts = mutate_all(epi_normal_raw_counts, function(x) as.numeric(as.character(x)))
epi_und_TME_raw_counts = mutate_all(epi_und_TME_raw_counts, function(x) as.numeric(as.character(x)))

#run full CNA pipeline - change counts and norm_cell depending on dataset
results <- pipelineCNA(epi_und_TME_raw_counts,
                    sample = 'tumor_class_w_normal',
                       SUBCLONES = TRUE,
                       ClonalCN = TRUE,
                       plotTree = TRUE
                       )

table(results$class)

#output to integrate in python
write.csv(results,"SCEVAN_results_epi_und_TME.csv", row.names = TRUE)

#load results with various cell type inclusions
results_epi_noRef = read.csv('SCEVAN_results_epi_noRef.csv', row.names='X')
results_epi_und_noRef = read.csv('SCEVAN_results_epi_und_noRef.csv', row.names='X')
results_epi_normal_ref = read.csv('SCEVAN_results_epi_normal_ref.csv', row.names='X')
results_epi_und_ref = read.csv('SCEVAN_results_epi_und_Ref.csv', row.names='X')

#check number of tumor cells classified by each
table(results_epi_noRef$class)
table(results_epi_und_noRef$class)
table(results_epi_normal_ref$class)
table(results_epi_und_ref$class)

#get id of all cells classified as tumor with each data subset
tum_ids_epi_only = rownames(results_epi_noRef[which(results_epi_noRef$class=='tumor'),])
tum_ids_epi_und = rownames(results_epi_und_noRef[which(results_epi_und_noRef$class=='tumor'),])
tum_ids_epi_normal = rownames(results_epi_normal_ref[which(results_epi_normal_ref$class=='tumor'),])
tum_ids_epi_und_normal = rownames(results_epi_und_ref[which(results_epi_und_ref$class=='tumor'),])

#check for overlap between these methods
length(tum_ids_epi_only[which(tum_ids_epi_only %in% tum_ids_epi_normal)])
length(tum_ids_epi_normal[which(tum_ids_epi_normal %in% tum_ids_epi_und)])
length(tum_ids_epi_und[which(tum_ids_epi_und %in% tum_ids_epi_und_normal)])
length(tum_ids_epi_und_normal[which(tum_ids_epi_und_normal %in% tum_ids_epi_normal)])

