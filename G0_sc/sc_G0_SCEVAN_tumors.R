#### Calculating G0 score for sc tumour data
library(readxl)
library(ggplot2)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

#retrieve classications as tumour/not from SCEVAN results
results_noRef = read.csv("SCEVAN_results.csv")

#retrieve classifications w/ healthy reference
results_ref = read.csv("SCEVAN_results_wREF.csv")


#need Cenk to run G0 scorer for me -------- tumor cells w/ reference
tumoural_cells <- epi_normal_counts[,colnames(epi_normal_counts) %in% results_ref$X[results_noRef$class=="tumor"]]
write.csv(tumoural_cells,"tumor_mat_wREF.csv")
##double check but looks like SCEVAN has identified some normal epithelial cells as tumors


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

names(sigs)[1] = 'G0sig.up'
names(sigs)[2] = 'G0sig.down'

#sigs[[2]] = sigs[[2]][sigs[[2]] %in% rownames(tumoural_cell_mat)]
#sigs[[1]] = sigs[[1]][sigs[[1]] %in% rownames(tumoural_cell_mat)] 

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

saveRDS(sigs, 'G0_sigs.RData')

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

#analysis classification - Quratile cutoffs -----------------


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

#analysis classification - more stringent cutoffs -----------------

slow_cut = -1
fast_cut = 1

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

write.csv(mat.scores, 'G0_scored_stringent.csv')

mat.scores

colnames(mat.scores)[7] = 'Cycle_type'

ggplot(data = mat.scores, aes(prolif_z, fill = Cycle_type)) +
  scale_color_manual(values=c("#bbbbbb", "#286837", "#772989"))+
  scale_fill_manual(values=c("#bbbbbb", "#286837", "#772989")) +
  geom_histogram(breaks = c(seq(-3.5,3.5, 0.25)),color="black") +
  theme_classic() +
  labs(
    x = "Proliferative capacity",
    y = "Counts"
  )


#checking signature dropout -----------------

upregulated = sigs$G0sig.up
downregulated = sigs$G0sig.down

up_tum_counts = tumoural_cells[which(row.names(tumoural_cells) %in% upregulated),]
down_tum_counts = tumoural_cells[which(row.names(tumoural_cells) %in% downregulated),]

sig_up0Props = rowMeans(up_tum_counts == 0)
sig_down0Props = rowMeans(down_tum_counts == 0)

sig_up0counts = rowSums(up_tum_counts == 0)
sig_down0counts = rowSums(down_tum_counts == 0)

hist(sig_up0Props)
hist(sig_down0Props)
barplot(sig_up0counts)
barplot(sig_down0counts)

## plot genes measured

length(upProps)
length(upProps[which(upProps > 250)])

cell_cut_off = 250

zero_Up_prop = length(up0counts[which(up0counts > cell_cut_off)])/length(up0counts)

zero_Down_prop = length(down0counts[which(down0counts > cell_cut_off)])/length(down0counts)

##### Dropout in genes not in the signature

non_sig_counts = tumoural_cells[-which(row.names(tumoural_cells) %in% downregulated),]
non_sig_counts = non_sig_counts[-which(row.names(non_sig_counts) %in% upregulated),]

nrow(tumoural_cells) - nrow(non_sig_counts)

nonSig_Props = rowMeans(non_sig_counts == 0)
nonSig_counts = rowSums(non_sig_counts == 0)


hist(nonSig_up0Props)
hist(nonSig_down0Props)
barplot(nonSig_up0Props)
barplot(nonSig_down0Props)

#### Comparing sig and not sig dropouts
# don't think this is really valid 
# - expect lots of genes to be downregulated/0 across the whole genome

mean(nonSig_up0Props)
mean(sig_up0Props)

mean(nonSig_down0Props)
mean(sig_down0Props)

df <- data.frame(
  type = c("sigUp", "sigDown", "nonSig"),
  prop   = c(mean(sig_up0Props), mean(sig_down0Props), mean(nonSig_Props))
)

df$type = as.factor(df$type)

ggplot(df, aes(x = type, y = prop, fill = type)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean gene dropout proportion in tumor cells", x = "Gene category", y = "Value") +
  theme_minimal()




