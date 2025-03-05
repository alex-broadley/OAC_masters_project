library(ggplot2)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

tum_categories = read.csv('G0_scored_tumorCells.csv',row.names=1)
tumor_counts = read.csv('tum_counts.csv', row.names = 1)

fast_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$strict_cycle_status == 'fast')])]
G0_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$strict_cycle_status == 'slow')])]
slow_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$strict_cycle_status == 'inter')])]

tumor_counts = t(tumor_counts)

tumor_counts = merge(as.data.frame(tumor_counts), as.data.frame(tum_categories)$strict_cycle_status, key = row)
tumor_counts$strict_cycle_status

tumor_counts = merge(tumor_counts, tum_categories, by = 0) 

tumor_counts[,'MKI67']
rownames(tumor_counts) = tumor_counts$Row.names

ggplot(tumor_counts, aes(x=strict_cycle_status, y=MKI67, fill = strict_cycle_status)) + 
  geom_violin() + 
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_minimal()

ggplot(tumor_counts, aes(x=strict_cycle_status, y=TOP2A, fill=strict_cycle_status)) + 
  geom_violin() + 
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_minimal()
