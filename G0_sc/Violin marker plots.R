library(ggplot2)
library(ggpubr)
library(rstatix)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

#tum_categories = read.csv('G0_scored_tumorCells.csv',row.names=1)
tum_categories = read.csv('G0_agreement_quartiles.csv',row.names=1)
epi_counts = read.csv("agreement_tum_SCRAN_counts.csv", row.names='X')

tumor_counts <- epi_counts[,colnames(epi_counts) %in% rownames(tum_categories)]

tumor_counts = t(tumor_counts)
tumor_counts = merge(tumor_counts, tum_categories, by = 0) 
table(tumor_counts$Cycle_type)

rownames(tumor_counts) = tumor_counts$Row.names

tumor_counts$Cycle_type = factor(tumor_counts$Cycle_type , levels=c("fast cycling", "cycling", "G0 arrested"))
tumor_counts$Cycle_type

tumor_counts$Cycle_type = str_to_sentence(tumor_counts$Cycle_type)

tumor_counts$Cycle_type = factor(tumor_counts$Cycle_type , levels=c("Fast cycling", "Cycling", "G0 arrested"))
tumor_counts$Cycle_type


my_comparisons = list(c("Fast cycling", "G0 arrested"), c("G0 arrested", 'Cycling'), c("Fast cycling", "Cycling"))

compare_means(MKI67 ~ Cycle_type,  data = tumor_counts)
compare_means(TOP2A ~ Cycle_type,  data = tumor_counts)

MKI67 = ggplot(tumor_counts, aes(x=Cycle_type, y=MKI67, fill = Cycle_type)) + 
  geom_violin() + 
  stat_compare_means(comparisons = my_comparisons, size = 5) +
  stat_compare_means(label.y = 4, size = 5) +
  scale_fill_manual(values=c("#e60000", "#b3b3cc", "#0000ff"), name = "Proliferation Type") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8, alpha = 0.3) +
  labs(x = '', title = 'MKI67', y = "Normalised MKI67 expression level") +
  theme_bw() +
  theme(axis.text=element_text(size=15), 
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none")

TOP2A = ggplot(tumor_counts, aes(x=Cycle_type, y=TOP2A, fill=Cycle_type)) + 
  geom_violin() + 
  stat_compare_means(comparisons = my_comparisons, size = 5) +
  stat_compare_means(label.y = 2.5, size = 5) +
  scale_fill_manual(values=c("#e60000", "#b3b3cc", "#0000ff"), name = "Proliferation Type") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8, alpha = 0.3) +
  labs(x = '', title = 'TOP2A', y = "Normalised TOP2A expression level") +
  theme_bw() +
  theme(axis.text=element_text(size=15), 
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

figure = grid.arrange(MKI67, TOP2A, ncol = 2, widths = c(0.8, 1))
annotate_figure(figure, top = text_grob("Comparison of proliferation markers in G0 arrested and fast cycling tumor cells", size = 20, face = "bold"))

