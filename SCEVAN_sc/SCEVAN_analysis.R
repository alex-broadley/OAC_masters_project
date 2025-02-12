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

#### Running G0 score:

#get expression matrix from only tumoural cells
tumoural_cells <- epi_counts[,colnames(epi_counts) %in% rownames(results_noRef[results_noRef$class=="tumor",])]
tumoural_cell_mat <- as.matrix(tumoural_cells)



setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')

library(readxl)

#read G0 arrest signature
G0sig <- read_xlsx("G0arrest_signature.xlsx")
G0sig.up <- G0sig[which(G0sig$Expression=="upregulated"),]$Gene
G0sig.down <- G0sig[which(G0sig$Expression=="downregulated"),]$Gene


load('downregulated_common.RData')
load('upregulated_common.RData')

sigs <- list(upregulated_common, downregulated_common)
#sigs[[1]] <- G0sig.up
#sigs[[2]] <- G0sig.down

#sigs[[2]] = sigs[[2]][sigs[[2]] %in% rownames(tumoural_cell_mat)]
#sigs[[1]] = sigs[[1]][sigs[[1]] %in% rownames(tumoural_cell_mat)] 

# Calculate G0 scores

library(GSVA)
## build GSVA parameter object
gsvapar <- gsvaParam(tumoural_cell_mat, sigs)


print(table(rownames(as.matrix(tumoural_cells)) %in% sigs[[1]]))
print(table(rownames(as.matrix(tumoural_cells)) %in% sigs[[2]]))


## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar, verbose = TRUE) #- doens't work for me for some reason

# load in results from Cenk
load('gsva_es.RData')

rownames(gsva_es) <- c("G0.up","G0.down")
mat.scores <- t(gsva_es)
mat.scores <- data.frame(mat.scores)
mat.scores$Sample <- rownames(mat.scores)
mat.scores$ProliferativeCapacity <- mat.scores$G0.down-mat.scores$G0.up
mat.scores$prolif_z = scale(mat.scores$ProliferativeCapacity, center = T, scale = T)

summary(mat.scores$prolif_z)
hist(mat.scores$prolif_z, breaks = seq(-3,3.5, 0.2))

mat.scores$cycle_status = 0

for (i in 1:nrow(mat.scores)){
  #paste(i, mat.scores$prolif_z)
  if (mat.scores$prolif_z[i] >= 0) {
    print('fast')
    mat.scores$cycle_status[i] = 'fast'
  } else {
    print("slow")
   mat.scores$cycle_status[i] = 'slow'
  }
}

mat.scores$cycle_status = factor(mat.scores$cycle_status, levels = c('slow', 'fast'))
table(mat.scores$cycle_status)
  
slow_cut = -0.75333
fast_cut = 0.64841



#couldn't get maxstat to generate two cutoffs, so using 1st and 3rd quartiles instead
library(maxstat)
maxstat.test(formula = cycle_status ~ prolif_z, data = mat.scores, pmethod="Lau92")

#Load RDAta from generated output file - don't think really need this 
setwd("/Users/alex/Documents/University/Year 4/BIOL0041/SC_analysis/output")
load("_count_mtx_annot.RData")
load("_CNAmtxSubclones.RData")
load("_CNAmtx.RData")
load(" _SubcloneDiffAnalysis.RData")


# Load the data
integrated <- readRDS(paste0(project_dir, "/02_G0arrestInMalignantCells/data/integrated_with_quiescence.rds"))
integrated <- DietSeurat(integrated, assays = "RNA")

# Subset the data to only include malignant cells
malignant <- subset(
  integrated,
  subset = type %in% c("TNBC", "ER", "PR", "HER", "DCIS", "neoplasm") &
    malignancy %in% "malignant"
)

rm(integrated)

# Run SCEVAN for ER G0 subtype-------------------------------------------------
er_g0 <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Quiescent" & type %in% "ER"
)

dir.create("ER_G0_arrested")
setwd("ER_G0_arrested")

# Get the counts
count_mtx <- GetAssayData(er_g0, layer = "counts", assay = "RNA")

set.seed(42)
if (ncol(count_mtx) > 3000) {
  count_mtx <- count_mtx[, sample(ncol(count_mtx), 3000)]
}

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "ER_G0_arrested",
  par_cores = 4,
  SUBCLONES = TRUE,
  plotTree = TRUE
) # this function generated all panels of Figure 3 and Extended Data Figure 3

setwd("..")

# Run SCEVAN for TNBC G0 subtype-----------------------------------------------
tnbc_g0 <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Quiescent" & type %in% "TNBC"
)

dir.create("TNBC_G0_arrested")
setwd("TNBC_G0_arrested")

# Get the counts
count_mtx <- GetAssayData(tnbc_g0, layer = "counts", assay = "RNA")

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "TNBC_G0_arrested",
  par_cores = 20,
  SUBCLONES = TRUE,
  plotTree = TRUE
)

setwd("..")

# Run SCEVAN for ER Fast cycling subtype---------------------------------------
er_fast <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Proliferating" & type %in% "ER"
)

dir.create("ER_Fast_cycling")
setwd("ER_Fast_cycling")

# Get the counts
count_mtx <- GetAssayData(er_fast, layer = "counts", assay = "RNA")

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "ER_Fast_cycling",
  par_cores = 20,
  SUBCLONES = TRUE,
  plotTree = TRUE
)

setwd("..")

# Run SCEVAN for TNBC Fast cycling subtype-------------------------------------
tnbc_fast <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Proliferating" & type %in% "TNBC"
)

dir.create("TNBC_Fast_cycling")
setwd("TNBC_Fast_cycling")

# Get the counts
count_mtx <- GetAssayData(tnbc_fast, layer = "counts", assay = "RNA")

set.seed(42)
if (ncol(count_mtx) > 3000) {
  count_mtx <- count_mtx[, sample(ncol(count_mtx), 3000)]
}

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "TNBC_Fast_cycling",
  par_cores = 20,
  SUBCLONES = TRUE,
  plotTree = TRUE
)

setwd("..")