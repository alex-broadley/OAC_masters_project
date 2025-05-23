## Feature importance plot

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

importances = read.csv('ordered_gene_importances.csv', row.names = 'X')

library(ggplot2)
library(tidyverse)

importances$Gene <- factor(importances$Gene, 
                              levels = importances$Gene[order(importances$mean_importance)])


ggplot(data = importances[1:25,], aes(x = Gene, y = mean_importance)) +
         geom_bar(stat = "identity", width = 0.65, fill = "#0000b3") + 
        theme_bw() +
        coord_flip() +
        labs(title = 'Mean importance per gene across RF and GB models', 
             x = 'Gene name',
             y = 'Mean importance') +
        theme(plot.title = element_text(hjust = 0.5))

print(importances[1:50, 'Gene'])
