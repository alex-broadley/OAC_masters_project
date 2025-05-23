library(DESeq2)
library(vsn)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(pathfindR)
library(patchwork)
library(EnhancedVolcano)
library(pathfindR)
library(ggpubr)

#### Fast vs G0 cell DGE

fast_G0_DE = read.csv('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data/fast_G0_DE.csv')

fast_G0_DE$padj <- as.numeric(fast_G0_DE$pvals_adj)
fast_G0_DE$log2FoldChange <- as.numeric(fast_G0_DE$logfoldchanges)


#run pathdindR w/ KEGG pathways
pathfindR_df = fast_G0_DE[, c('names', 'logfoldchanges', 'pvals_adj')]

sig_genes = pathfindR_df[which(pathfindR_df$pvals_adj < 0.05),]


pathfindR_df = pathfindR_df[which(pathfindR_df$pvals_adj < 0.05),]


#run pathdindR w/ KEGG pathways
output_df <- run_pathfindR(pathfindR_df,
                           output_dir="Fast_GO_DGE", gene_sets="KEGG")

od_SC <- data.frame(output_df)


ggplot(data = fast_G0_DE, aes(pvals_adj)) +
  theme_classic() +
  geom_histogram(fill="#0000b3", color='white') +
  labs(
    x = "adjust P-value",
    y = "Frequency",
    title = 'P-value histogram for DGE analysis comparing aggressive vs slowly proliferating primary tumors'
  ) +   theme(axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=15),
              axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15),
              plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
              legend.text = element_text(size = 10))


### Output DEG's as space separated for Gprofiler

writeLines(as.character(sig_genes$names), "SC_DGE_names.txt")





#### Volcano Plot ------------------------------------------------------------------------------

fast_G0_DE$sig_direction = "Not significant"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges > 0 & fast_G0_DE$pvals_adj < 0.05] <- "Significant"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges < 0 & fast_G0_DE$pvals_adj < 0.05] <- "Significant"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges > 0.5 & fast_G0_DE$pvals_adj < 0.05] <- "Significant & log2FC > 0.05"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges < -0.5 & fast_G0_DE$pvals_adj < 0.05] <- "Significant & log2FC < -0.05"

fast_G0_DE$delabel = NA

log2FC_high_df = fast_G0_DE[abs(fast_G0_DE$logfoldchanges) > 0.5 & fast_G0_DE$pvals_adj < 0.05, ]

for (row in 1:nrow(fast_G0_DE)) {
  if (fast_G0_DE$names[row] %in% head(log2FC_high_df[order(log2FC_high_df$pvals_adj), "names"], 50)){
    fast_G0_DE$delabel[row] = fast_G0_DE$names[row]
  }
  else if (fast_G0_DE$sig_direction[row] == 'Not significant'){
    fast_G0_DE$delabel[row] = NA
  }
  else if (fast_G0_DE$sig_direction[row] == 'Significant & log2FC < 0'){
    fast_G0_DE$delabel[row] = NA
  }
  else if (fast_G0_DE$sig_direction[row] == 'Significant & log2FC > 0'){
    fast_G0_DE$delabel[row] = NA
  }
  else {
    fast_G0_DE$delabel[row] = NA
  }
}


#### create volcano plot 
sc_volcano = ggplot(data = fast_G0_DE, aes(x = logfoldchanges, y = -log10(pvals_adj), colour = sig_direction, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  geom_point() +
  scale_color_manual(values = c("lightgrey", "darkgrey", "darkblue", "red"), name = "") +
  labs(x = "log2FC", y = "-log10(p-adj)", title = "DEG Volcano Plot") +
  theme_classic() +
  geom_text_repel(max.overlaps = 12) + 
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
        legend.text = element_text(size = 10))

sc_volcano = ggpar(sc_volcano, legend = "top")

## Create enrichment bar plot

sc_enrich_bar = ggplot(data = od_SC[0:10,], aes(x = Term_Description, y = Fold_Enrichment, fill = highest_p)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(x = "", y = "Fold Enrichment", title = "KEGG pathway enrichment (top 10 pathways)") +
  theme_classic() +
  scale_fill_gradient(low = '#b3908b', high = 'darkred', (name = 'P-value')) + 
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=0),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10))

sc_enrich_bar = ggpar(sc_enrich_bar, legend = "right")


library(gridExtra)

figure = grid.arrange(sc_volcano, sc_enrich_bar, ncol = 2, widths = c(1, 1.1))
annotate_figure(figure, top = text_grob("Summary of DGE results in fast cylcing vs G0 arrested tumor cells", size = 20, face = "bold"))


diff_selected_genes = c('TP53', 'CDKN2A', 'SMAD4', 'SPG20', 'FBXW7', 'KRAS', 'ARID1A', 'COL19A1', 'DHX16')


sig_genes[sig_genes$names %in% diff_selected_genes,]


