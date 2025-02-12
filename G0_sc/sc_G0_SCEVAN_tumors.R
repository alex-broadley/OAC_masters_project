#### Calculating G0 score for sc tumour data
library(readxl)

#retrieve classications as tumour/not from SCEVAN results
results_noRef = read.csv("SCEVAN_results.csv")

#get expression matrix from only tumoural cells
tumoural_cells <- epi_counts[,colnames(epi_counts) %in% rownames(results_noRef[results_noRef$class=="tumor",])]
tumoural_cell_mat <- as.matrix(tumoural_cells)

setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')

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