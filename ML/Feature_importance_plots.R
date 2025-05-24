### Feature importance plots

library(ggplot2)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

#read in bulk and single cell averaged ML feature importances
bulk_importances = read.csv('Bulk_feature_importances.csv', row.names = 'X')
SC_importances = read.csv('SC_feature_importance.csv', row.names = 'X')

bulk_importances$gene_name <- factor(bulk_importances$gene_name, levels = bulk_importances$gene_name[order(bulk_importances$mean_importance, decreasing = FALSE)])
SC_importances$gene_name <- factor(SC_importances$gene_name, levels = SC_importances$gene_name[order(SC_importances$mean_importance, decreasing = FALSE)])


#create stick plot of feature importances - BULK
bulk_feature_importances = ggplot(bulk_importances[0:15,], aes(x = mean_importance, y = gene_name)) +
  geom_point(color = '#3366ff') +
  labs(
    x = "Mean feature importance",
    y = "Gene",
    title = "RF-XGB averaged feature importances - Bulk RNA"
  ) +
  theme_classic() +
  geom_segment(aes(x = 0, xend = mean_importance, y = gene_name, yend = gene_name)) +
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5, face = 'bold'))

#create stick plot of feature importances - single cell
SC_feature_importances = ggplot(SC_importances[0:15,], 
                                  aes(x = mean_importance, y = gene_name)) +
  geom_point(color = '#3366ff') +
  labs(
    x = "Mean feature importance",
    y = "Gene",
    title = "RF-XGB averaged feature importances - Single cell"
  ) +
  theme_classic() +
  geom_segment(aes(x = 0, xend = mean_importance, y = gene_name, yend = gene_name)) +
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5, face = 'bold'))


### Elbow plots to determine cutoffs for functional enrichment analysis -----------------------------------

ggplot(data = bulk_importances, aes(mean_importance)) +
  geom_histogram(breaks = c(seq(0,0.02, 0.0005)), color="white", fill = "darkblue") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Proliferative capacity",
    y = "Counts",
    title = "Bulk feature importances"
  ) +
  geom_vline(xintercept = 0.0025, col = "red", linetype = 'dashed') +
  theme(axis.text=element_text(size=15), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  labs(
    x = "Mean feature importances",
    y = "Counts",
    title = "Primary tumor mean feature importance distribution",
  )

ggplot(data = SC_importances, aes(mean_importance)) +
  geom_histogram(breaks = c(seq(0,0.02, 0.0005)), color="white", fill = "darkblue") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0.0035, col = "red", linetype = 'dashed') +
  theme(axis.text=element_text(size=15), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  labs(
    x = "Mean feature importances",
    y = "Counts",
    title = "Single Cell mean feature importance distribution",
  )

### Extracting gene lists based on elbow -----------------------------------

bulk_important_genes = as.character(list(bulk_importances$gene_name[bulk_importances$mean_importance > 0.004])[[1]])
SC_important_genes = as.character(list(SC_importances$gene_name[SC_importances$mean_importance > 0.0035])[[1]])

cat(bulk_important_genes, file = "bulk_ML_important_genes.txt", sep = " ")
cat(SC_important_genes, file = "SC_ML_important_genes.txt", sep = " ")

### Extracting top 20 most important-------------------------------------------------------

top_20_bulk = bulk_importances[order(bulk_importances$mean_importance, decreasing = TRUE), ][1:20, ]
top_20_sc = SC_importances[order(SC_importances$mean_importance, decreasing = TRUE), ][1:20, ]


write.table(top_20_bulk$gene_name, file = "top_20_bulk.txt", sep = " ", row.names = FALSE, quote = FALSE)
write.table(top_20_sc$gene_name, file = "top_20_sc", sep = " ", row.names = FALSE, quote = FALSE)




