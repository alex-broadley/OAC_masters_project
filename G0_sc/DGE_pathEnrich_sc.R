#Initialise packages
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(pathfindR)
library(patchwork)
library(pathfindR)
library(ggpubr)

#### Fast vs G0 cell DGE

#read scanpy 'rank_genes_groups' DEGs between fast and G0 arrested tumor cells
fast_G0_DE = read.csv('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data/fast_G0_DE.csv')

#ensure pvals and logfold changes are numeric
fast_G0_DE$padj <- as.numeric(fast_G0_DE$pvals_adj)
fast_G0_DE$log2FoldChange <- as.numeric(fast_G0_DE$logfoldchanges)

#get relevant columns for DEGs (p-adj <0.05)
pathfindR_df = fast_G0_DE[, c('names', 'logfoldchanges', 'pvals_adj')]
pathfindR_df = pathfindR_df[which(pathfindR_df$pvals_adj < 0.05),]

#run pathdindR w/ KEGG pathways
output_df <- run_pathfindR(pathfindR_df,
                           output_dir="Fast_GO_DGE", gene_sets="KEGG")
od_SC <- data.frame(output_df)

#DEPRECATED - p-value hist checking for 'anti-conservative', but do not expect this distribution in single cell DGEs
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

#### Volcano Plot ------------------------------------------------------------------------------

# Create labels for DEGs based on signficance and direction of expression
fast_G0_DE$sig_direction = "Not significant"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges > 0 & fast_G0_DE$pvals_adj < 0.05] <- "Significant"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges < 0 & fast_G0_DE$pvals_adj < 0.05] <- "Significant"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges > 0.5 & fast_G0_DE$pvals_adj < 0.05] <- "Significant & log2FC > 0.05"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges < -0.5 & fast_G0_DE$pvals_adj < 0.05] <- "Significant & log2FC < -0.05"


#set which genes want to label on plot - here choosing 50 with highest pvals and abs(LFC) > 0.5
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


#### create volcano plot with same theme choices as in all other plots
sc_volcano = ggplot(data = fast_G0_DE, aes(x = logfoldchanges, y = -log10(pvals_adj), colour = sig_direction, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  geom_point() +
  scale_color_manual(values = c("lightgrey", "darkgrey", "blue", "red"), name = "") +
  labs(x = "log2FC", y = "-log10(p-adj)", title = "DEG Volcano Plot") +
  theme_classic() +
  geom_text_repel(max.overlaps = 12,
                  min.segment.length = unit(0, 'lines')) + 
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
        legend.text = element_text(size = 10))

#move legend to top
sc_volcano = ggpar(sc_volcano, legend = "top")


### Bar plot for enrichment results
#plot 10 most significant pathways followed, color by highest p-value
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

#move legend to top of figure
sc_enrich_bar = ggpar(sc_enrich_bar, legend = "right")



#put volcano and enrichment bar chart next to each other
library(gridExtra)
figure = grid.arrange(sc_volcano, sc_enrich_bar, ncol = 2, widths = c(1, 1.1))
#annotate figure with title
annotate_figure(figure, top = text_grob("Summary of DGE results in fast cylcing vs G0 arrested tumour cells", size = 20, face = "bold"))




