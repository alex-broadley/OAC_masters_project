#### Calculating G0 score for sc tumour data
library(readxl)
library(ggplot2)
library(tidyverse)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

#retrieve classications as tumour/not from SCEVAN results ------------------------

#results_epi_noRef = read.csv('SCEVAN_results_epi_noRef.csv', row.names='X')
#results_epi_und_noRef = read.csv('SCEVAN_results_epi_und_noRef.csv', row.names='X')
#results_epi_normal_ref = read.csv('SCEVAN_results_epi_normal_ref.csv', row.names='X')
#results_epi_und_ref = read.csv('SCEVAN_results_epi_und_Ref.csv', row.names='X')


#SCEVAN_results_toUse = results_epi_und_noRef

#discuss w/ Maria but seems to be identifying some of the 'normal' cells as tumor cells
#tumor_scran_counts = read.csv("epi_und_SCRAN_counts.csv", row.names='X')
tumor_scran_counts = read.csv('agreement_tum_SCRAN_counts.csv', row.names='X')

#tumor_scran_counts <- tumor_scran_counts[,colnames(tumor_scran_counts) %in% rownames(SCEVAN_results_toUse)[SCEVAN_results_toUse$class=="tumor"]]
tumor_scran_mat <- as.matrix(tumor_scran_counts)

#check this is 912 cells as expected
length(colnames(tumor_scran_mat))

#check these are SCRAN normalised counts
max(tumor_scran_mat)

#read G0 arrest signature ------------------------------------------------------
setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')


load('downregulated_common.RData')
load('upregulated_common.RData')

sigs <- list(upregulated_common, downregulated_common)
names(sigs)[1] = 'G0sig.up'
names(sigs)[2] = 'G0sig.down'


# Calculate G0 scores ----------------------------------------------------------

# only retain variable genes
keep <- apply(tumor_scran_mat, 1, function(x) var(x, na.rm = T) > 0)
tumor_scran_mat <- data.matrix(tumor_scran_mat[keep, ])



library(GSVA)
## build GSVA parameter object
gsvapar <- gsvaParam(tumor_scran_mat, sigs)

print(table(rownames(as.matrix(tumor_scran_mat)) %in% sigs[[1]]))
print(table(rownames(as.matrix(tumor_scran_mat)) %in% sigs[[2]]))

sigs[[2]] = sigs[[2]][sigs[[2]] %in% rownames(as.matrix(tumor_scran_mat))]

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar, verbose = TRUE) 

rownames(gsva_es)

mat.scores = t(gsva_es)
mat.scores <- data.frame(mat.scores)
mat.scores$ProliferativeCapacity <- mat.scores$G0sig.down - mat.scores$G0sig.up
mat.scores$prolif_z = scale(mat.scores$ProliferativeCapacity, center = T, scale = T)

summary(mat.scores$prolif_z)
hist(mat.scores$prolif_z)

mat.scores$cycle_status = 0

#split on 0 ----------------------------------------------------

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

mat.scores = read.csv('G0_agreement_quartiles.csv')

#analysis classification - Quratile cutoffs -----------------------------------
slow_cut = quantile(mat.scores$prolif_z, probs = c(0.25,0.75))[[1]]
fast_cut = quantile(mat.scores$prolif_z, probs = c(0.25,0.75))[[2]]

for (i in 1:nrow(mat.scores)){
  if (mat.scores$prolif_z[i] >= fast_cut) {
    mat.scores$Cycle_type[i] = 'Fast cycling'
  } else if (mat.scores$prolif_z[i] <= slow_cut) {
    mat.scores$Cycle_type[i] = 'G0 arrested'
  } else {
    mat.scores$Cycle_type[i] = 'Cycling'
  }
}

table(mat.scores$Cycle_type)

write.csv(mat.scores, 'G0_agreement_quartiles.csv')

mat.scores

mean(mat.scores$prolif_z)

ggplot(data = mat.scores, aes(prolif_z, fill = Cycle_type)) +
  scale_color_manual(values=c("#bbbbbb", "#e60039", "#3366ff"))+
  scale_fill_manual(values=c("#bbbbbb", "#e60039", "#3366ff")) +
  geom_histogram(breaks = c(seq(-15,7, 0.25), slow_cut, fast_cut),color="black") +
  theme_classic() +
  theme(axis.text=element_text(size=15), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  labs(
    x = "Proliferative score",
    y = "Counts",
    title = "Tumor cell proliferative score distribution",
    fill = 'Proliferation type'
  )

#analysis classification - more stringent cutoffs -----------------

slow_cut = -1
fast_cut = 1

# try with 1 and -1 cutoffs

for (i in 1:nrow(mat.scores)){
  #paste(i, mat.scores$prolif_z)
  if (mat.scores$prolif_z[i] >= fast_cut) {
    mat.scores$Cycle_type[i] = 'fast cycling'
  } else if (mat.scores$prolif_z[i] <= slow_cut) {
    mat.scores$Cycle_type[i] = 'G0 arrested'
  } else {
    mat.scores$Cycle_type[i] = 'cycling'
  }
}

table(mat.scores$Cycle_type)

write.csv(mat.scores, 'G0_agreement_stringent.csv', row.names = TRUE)

mat.scores

ggplot(data = mat.scores, aes(prolif_z, fill = Cycle_type)) +
  scale_color_manual(values=c("#bbbbbb", "#286837", "#772989"))+
  scale_fill_manual(values=c("#bbbbbb", "#286837", "#772989")) +
  geom_histogram(breaks = c(seq(-15,7, 0.25)),color="black") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Proliferative capacity",
    y = "Counts",
    title = "Cycle type with stringent cutoffs"
  )


#checking signature dropout -----------------

upregulated = sigs$G0sig.up
downregulated = sigs$G0sig.down

tumoural_cells = tumor_scran_mat

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

cell_cut_off = 250

zero_Up_prop = length(sig_up0counts[which(sig_up0counts > cell_cut_off)])/length(sig_up0counts)
zero_Down_prop = length(sig_down0counts[which(sig_down0counts > cell_cut_off)])/length(sig_down0counts)

##### Dropout in genes not in the signature

non_sig_counts = tumoural_cells[-which(row.names(tumoural_cells) %in% downregulated),]
non_sig_counts = non_sig_counts[-which(row.names(non_sig_counts) %in% upregulated),]

nrow(tumoural_cells) - nrow(non_sig_counts)

nonSig_Props = rowMeans(non_sig_counts == 0)
nonSig_counts = rowSums(non_sig_counts == 0)


hist(nonSig_Props)
hist(nonSig_counts)
barplot(nonSig_Props)
barplot(nonSig_counts)

#### Comparing sig and not sig dropouts
# don't think this is really valid 
# - expect lots of genes to be downregulated/0 across the whole genome
# could remove genes with 100% dropout rate -> not expressed in any cells?

mean(nonSig_Props)
mean(sig_up0Props)

mean(sig_up0counts)
mean(sig_down0counts)
mean(nonSig_counts)

df <- data.frame(
  type = c("sigUp", "sigDown", "nonSig"),
  prop   = c(mean(sig_up0Props), mean(sig_down0Props), mean(nonSig_Props))
)

df$type = as.factor(df$type)

ggplot(df, aes(x = type, y = prop, fill = type)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean gene dropout proportion in tumor cells", x = "Gene category", y = "Value") +
  theme_minimal()


## mean and sd expression level plot
library(dplyr)

# Calculate mean and standard deviation for each gene in tumor mat scores
mean_sd_counts <- tumor_scran_counts %>%
  rowwise() %>%
  mutate(
    Mean = mean(c_across(everything())),
    SD = sd(c_across(everything()))
  ) %>%
  ungroup() %>%
  select(Mean, SD)

mean_sd_counts <- as.data.frame(mean_sd_counts[, c("Mean", "SD")])
row.names(mean_sd_counts) <- row.names(tumor_scran_counts)


sigs

mean_sd_counts$signature = "No signature"
mean_sd_counts$signature[which(row.names(mean_sd_counts) %in% sigs$G0sig.up)] = 'Upregulated'
mean_sd_counts$signature[which(row.names(mean_sd_counts) %in% sigs$G0sig.down)] = 'Downregulated'

mean_sd_counts$signature = factor(mean_sd_counts$signature, levels = c('Downregulated', 'Upregulated', 'No signature'))

ggplot(mean_sd_counts[mean_sd_counts$signature != 'No signature',], aes(x=SD, y=Mean, group = signature)) +
  geom_point(size=1, alpha = 0.8, aes(color = signature)) +
  scale_color_manual(values=c('Red','Green', 'Black')) +
  theme_bw()


library(maxstat)
?maxstat
  