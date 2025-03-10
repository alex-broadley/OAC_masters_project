#### Calculating G0 score for sc tumour data
library(readxl)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

#retrieve classications as tumour/not from SCEVAN results
results_noRef = read.csv("SCEVAN_results.csv")
#read in epithelial cell caount matrix
epi_counts = read.csv("epi_raw_counts_matrix.csv")
rownames(epi_counts) = epi_counts$X
epi_counts <- subset(epi_counts, select = -X)

#get expression matrix from only tumoural cells
tumoural_cells <- epi_counts[,colnames(epi_counts) %in% results_noRef$X[results_noRef$class=="tumor"]]

tumoural_cell_mat <- as.matrix(tumoural_cells)
colnames(tumoural_cell_mat)
rownames(tumoural_cell_mat)

#write tumor expression matrix
write.csv(tumoural_cells,"tumor_expression_mat.csv")

setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')

#read G0 arrest signature
load('downregulated_common.RData')
load('upregulated_common.RData')

sigs <- list(upregulated_common, downregulated_common)
sigs[[1]] <- G0sig.up
sigs[[2]] <- G0sig.down

#sigs[[2]] = sigs[[2]][sigs[[2]] %in% rownames(tumoural_cell_mat)]
#sigs[[1]] = sigs[[1]][sigs[[1]] %in% rownames(tumoural_cell_mat)] 

saveRDS(sigs, 'G0_sigs.RData')

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

# Calculate G0 scores

library(GSVA)
## build GSVA parameter object
gsvapar <- gsvaParam(tumoural_cell_mat, sigs)

print(table(rownames(as.matrix(tumoural_cells)) %in% sigs[[1]]))
print(table(rownames(as.matrix(tumoural_cells)) %in% sigs[[2]]))

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar, verbose = TRUE) #- doesn't work for me for some reason

# load in results from Cenk
load('mat.scores.RData')

colnames(mat.scores) <- c("G0.up","G0.down", "Sample")
#mat.scores <- t(mat.scores)
mat.scores <- data.frame(mat.scores)
mat.scores$ProliferativeCapacity <- mat.scores$G0.down - mat.scores$G0.up
mat.scores$prolif_z = scale(mat.scores$ProliferativeCapacity, center = T, scale = T)

summary(mat.scores$prolif_z)
hist(mat.scores$prolif_z, breaks = seq(-3,3.5, 0.2))

mat.scores$cycle_status = 0

#exploratory classification

for (i in 1:nrow(mat.scores)){
  #paste(i, mat.scores$prolif_z)
  if (mat.scores$prolif_z[i] >= 0) {
    mat.scores$cycle_status[i] = 'fast'
  } else {
    mat.scores$cycle_status[i] = 'slow'
  }
}

mat.scores$cycle_status = factor(mat.scores$cycle_status, levels = c('slow', 'fast'))
table(mat.scores$cycle_status)

#analysis classification

slow_cut = -0.753605
fast_cut = 0.683496

# try with 1 and -1 cutoffs

for (i in 1:nrow(mat.scores)){
  #paste(i, mat.scores$prolif_z)
  if (mat.scores$prolif_z[i] >= fast_cut) {
    mat.scores$strict_cycle_status[i] = 'fast cycling'
  } else if (mat.scores$prolif_z[i] <= slow_cut) {
    mat.scores$strict_cycle_status[i] = 'G0 arrested'
  } else {
    mat.scores$strict_cycle_status[i] = 'cycling'
  }
}

table(mat.scores$strict_cycle_status)

write.csv(mat.scores, 'G0_scored_tumorCells.csv')

mat.scores

colnames(mat.scores)[7] = 'Cycle_type'

ggplot(data = mat.scores, aes(prolif_z, fill = Cycle_type)) +
  scale_color_manual(values=c("#bbbbbb", "#286837", "#772989"))+
  scale_fill_manual(values=c("#bbbbbb", "#286837", "#772989")) +
  geom_histogram(breaks = c(seq(-3.5,3.5, 0.25), slow_cut, fast_cut),color="black") +
  theme_classic() +
  labs(
    x = "Proliferative capacity",
    y = "Counts"
  )




