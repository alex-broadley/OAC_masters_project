library(ggplot2)
library(ggpubr)
library(rstatix)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

#tum_categories = read.csv('G0_scored_tumorCells.csv',row.names=1)
tum_categories = read.csv('G0_scored_stringent.csv',row.names=1)
tumor_counts = read.csv('tum_counts.csv', row.names = 1)

#fast_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$Cycle_type == 'fast')])]
#G0_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$Cycle_type == 'slow')])]
#slow_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$Cycle_type == 'inter')])]
names(tum_categories)[7] = "Cycle_type"

tumor_counts = t(tumor_counts)
tumor_counts = merge(tumor_counts, tum_categories, by = 0) 
table(tumor_counts$Cycle_type)

#tumor_counts$Cycle_type = gsub('slow', 'G0_arrested',tumor_counts$Cycle_type)
#tumor_counts$Cycle_type = gsub('inter', 'slow',tumor_counts$Cycle_type)


rownames(tumor_counts) = tumor_counts$Row.names

tumor_counts$Cycle_type = factor(tumor_counts$Cycle_type , levels=c("fast cycling", "cycling", "G0 arrested"))
tumor_counts$Cycle_type

my_comparisons = list(c("fast cycling", "G0 arrested"), c("G0 arrested", 'cycling'), c("fast cycling", "cycling"))

stat.test <- tumor_counts %>% 
  wilcox_test(MKI67 ~ Cycle_type, comparisons = my_comparisons) %>%
  add_xy_position(x = "Cycle_type") %>%
  mutate(myformatted.p = paste0("p = ", p))
stat.test

compare_means(MKI67 ~ Cycle_type,  data = tumor_counts)
compare_means(TOP2A ~ Cycle_type,  data = tumor_counts)

ggplot(tumor_counts, aes(x=Cycle_type, y=MKI67, fill = Cycle_type)) + 
  geom_violin() + 
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 30) +
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_minimal()

ggplot(tumor_counts, aes(x=Cycle_type, y=TOP2A, fill=Cycle_type)) + 
  geom_violin() + 
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 10) +
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_minimal()

