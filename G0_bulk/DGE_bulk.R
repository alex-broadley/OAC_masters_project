## function initialisation
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
library(pathfindR)




# creates a plot displaying the standard deviation across - adapated from AdCompBio practical
create_SDplot <- function(data, plot_title){
  plot = meanSdPlot(data, ranks = FALSE, plot = F)
  plot$gg +
    theme_bw() +
    labs(title = plot_title) +
    xlab('Mean read count')+
    ylab('Standard deviation') +
    theme(
      legend.position = "right",
      legend.key.size = unit(2, "cm"), 
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      axis.text.x = element_text(size = 25),  
      axis.text.y = element_text(size = 25),
      axis.title=element_text(size=30),
      plot.title = element_text(size = 40)
    ) 
}


## pre-processing --------------------------------------------------------------------

# load bulk RNA data files
setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')

#load sample ID labels
load("oac.expr.final.RData")

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

# load bulk file w/ attached G0 classifications -> take both the raw counts file and the normalised as classifications in normalised
bulk_RNA = read.csv('bulk_RNA_classified.csv', row.names = 'X')
bulk_RNA_rawCounts = read.delim('expressionMatrix.OCCAMS.RawCounts.txt')

# map the G0 classifications to the DNAid
classifications_mapped = bulk_RNA[,c("DNAid", "proliferation")]

# Make sample IDs a column for easier merging
bulk_RNA_rawCounts$Sample = rownames(bulk_RNA_rawCounts)

# create mapping between the DNAid and sample id
mapping <- unique(oac.expr[,c("Sample","DNAid")])

# map the DNA id to the sample ID in the raw counts file
mat.scores.mapped <- merge(bulk_RNA_rawCounts, mapping, by.x='Sample', by.y="Sample", all.x=F, all.y=F)
# map the DNAid to the classification in the raw counts file
raw_counts_mapped = merge(mat.scores.mapped, classifications_mapped, by.X = 'DNAid', by.y = 'DNAid', all.x = F, all.y = F)

# remove any non-classified tumors 
raw_counts_mapped = na.omit(raw_counts_mapped)
rownames(raw_counts_mapped) = raw_counts_mapped$DNAid

rm_col <- c("Sample", "DNAid", "proliferation")
raw_counts_final = raw_counts_mapped[, !(colnames(raw_counts_mapped) %in% rm_col), drop = FALSE]

raw_counts_final <- sapply(raw_counts_final, as.numeric)
rownames(raw_counts_final) = rownames(raw_counts_mapped)

# extract metadata for DGE
metadata = raw_counts_mapped[, c('DNAid', 'proliferation')]

# set factors for metadata - make output changes in fast relative to slow tumors
metadata$proliferation <- factor(metadata$proliferation)
metadata$proliferation <- relevel(metadata$proliferation, "slow")
metadata$DNAid <- factor(metadata$DNAid)

## DGE --------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = round(t(raw_counts_final)), colData = metadata, 
                              design = ~ proliferation)

#remove genes with read count <= 10
dds <- dds[rowSums(counts(dds)) > 10, ]

## QC ---------------------------------------------
#calculate size factors
dds <- estimateSizeFactors(dds)


#VSD transfom data
vsd <- vst(dds, blind = FALSE)

#generate plot for each performed transformation
SdPlot_vst  <- create_SDplot(assay(vsd), "VST transformation")

#get distance matrix using VST tranformed data
sampleDists_vsd <- dist(t(assay(vsd)))

#convert to matrix datatype
sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)

#create color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

#run PCA and create plot
PCA_plot <- plotPCA(vsd, intgroup=c("proliferation"))  +
  theme_bw()  +
  scale_color_brewer(palette = "Dark2")

PCA_plot 

### Running DESeq DGE pipeline --------------------------------------------------

dds_DGE <- DESeq(dds)

#extract results
dds_DGE_results <- results(dds_DGE)

#p-value histogram to check for 'anti-conservative distribution'
ggplot(data = dds_DGE_results, aes(padj)) +
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

#subset only the significant DEGs
sig_results = subset(dds_DGE_results, padj < 0.05)



sig_results[rownames(sig_results) == 'RPS2', ]


#get names of 50 genes with lowest P-values
selected_genes <- rownames(sig_results[order(sig_results$padj), ])[1:50]
#use name list to get expression levels for these genes
selected_expression <- assay(vsd)[selected_genes, ]
#scale expression for better heatmap plot
selected_expression <- selected_expression - rowMeans(selected_expression)

#outputting significant DEGs

df <- data.frame(sig_results[,
                             c("log2FoldChange","padj")])
df$Gene = rownames(df)

write.csv(df, "slow_v_fast_BULK_DEG_results.csv")


#### Pathway enrichment on results -----------------------------------------------

#make padj and log2fold columns numeric
df$padj <- as.numeric(df$padj)
df$log2FoldChange <- as.numeric(df$log2FoldChange)

#run pathdindR w/ KEGG pathways
output_df <- run_pathfindR(df[,c("Gene","log2FoldChange","padj")],
                           output_dir="slow_vs_fast_BULK", gene_sets="KEGG")
od <- data.frame(output_df)

#create term gene graph from the results
term_gene_graph(result_df = od, use_description = TRUE)

### Volcano Plot -----------------------------------------------------------------

volcano_df = dds_DGE_results
volcano_df$Gene = rownames(volcano_df)

# Create labels for DEGs based on signficance and direction of expression
volcano_df$sig_direction = "Not significant"
volcano_df$sig_direction[volcano_df$log2FoldChange > 0 & volcano_df$padj < 0.05] <- "Significant"
volcano_df$sig_direction[volcano_df$log2FoldChange < 0 & volcano_df$padj < 0.05] <- "Significant"
volcano_df$sig_direction[volcano_df$log2FoldChange > 0.5 & volcano_df$padj < 0.05] <- "Significant & log2FC > 0.05"
volcano_df$sig_direction[volcano_df$log2FoldChange < -0.5 & volcano_df$padj < 0.05] <- "Significant & log2FC < -0.05"

table(volcano_df$sig_direction)


#create df of genes with highest L2FC and significant
log2FC_high_df = volcano_df[abs(volcano_df$log2FoldChange) > 0.5, ]
log2FC_high_df <- na.omit(log2FC_high_df)
log2FC_high_df = log2FC_high_df[log2FC_high_df$padj < 0.05, ] 


#set which genes want to label on plot - here choosing 20 with highest pvals and abs(LFC) > 0.5
#also include label if one of the DGE's highlight in discussion
volcano_df$DEGlabel = NA

for (row in 1:nrow(volcano_df)) {
  if (volcano_df$Gene[row] %in% head(log2FC_high_df[order(log2FC_high_df$padj), "Gene"], 20)){
    volcano_df$DEGlabel[row] = volcano_df$Gene[row]
  }
  else if (volcano_df$sig_direction[row] == 'Not significant'){
    volcano_df$DEGlabel[row] = NA
  }
  else if (volcano_df$sig_direction[row] == 'Significant & log2FC < 0'){
    volcano_df$DEGlabel[row] = NA
  }
  else if (volcano_df$sig_direction[row] == 'Significant & log2FC > 0'){
    volcano_df$DEGlabel[row] = NA
  }
  else {
    volcano_df$DEGlabel[row] = NA
  }
}

#find indexes of genes want to highlight for discussion
which(volcano_df$Gene %in% c('SPARCL1'))
which(volcano_df$Gene %in% c('RPS2'))
which(volcano_df$Gene %in% c('FBXW7'))
which(volcano_df$Gene %in% c('TP53'))
which(volcano_df$Gene %in% c('TPX2'))
which(volcano_df$Gene %in% c('COL19A1'))

#manually assing DEG labels for now
volcano_df$DEGlabel[39411] = 'SPARCL1'
volcano_df$DEGlabel[36951] = 'RPS2'
volcano_df$DEGlabel[19833] = 'FBXW7'
volcano_df$DEGlabel[41357] = 'TP53'
volcano_df$DEGlabel[41458] = 'TPX3'
volcano_df$DEGlabel[16810] = 'COL19A1'


#### create volcano plot 
bulk_volcano = ggplot(data = volcano_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig_direction, label = DEGlabel)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(max.overlaps = 15, force = 10, colour = 'black',
                  min.segment.length = unit(0, 'lines')) + 
  labs(x = "log2FC", y = "-log10(p-adj)", title = "DEG Volcano Plot") +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  scale_color_manual(values = c("grey", "darkgrey", "blue", "red"), name = "") +
  theme_classic() +
  theme(axis.text.y = element_text(size=12),
                                           axis.text.x = element_text(size=15),
                                           axis.title.x = element_text(size=15),
                                           axis.title.y = element_text(size=15),
                                           plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
                                           legend.text = element_text(size = 10))

library(ggpubr)
bulk_volcano = ggpar(bulk_volcano, legend = "top")

#### Checking for DGE of differentially selected genes -------------------------------------------------------

diff_selected_genes = c('TP53', 'CDKN2A.p16INK4a', 'SMAD4', 'SPG20', 'FBXW7', 'KRAS', 'ARID1A', 'COL19A1', 'DHX16')

sig_results[rownames(sig_results) %in% diff_selected_genes,]

sig_results[rownames(sig_results) %in% 'RPS2',]


### Creating sideways bar plot for pathway enrichment ---------------------------------------------------------

bulk_enrich_bar = ggplot(data = od[0:10,], aes(x = Term_Description, y = Fold_Enrichment, fill = highest_p)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(x = "", y = "Fold Enrichment", title = "KEGG pathway enrichment (top 10 pathways)") +
  theme_classic() +
  scale_fill_gradient(low = '#b3908b', high = 'darkred', (name = 'P-adj')) + 
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=0),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 7))

#put legend at top of plot
bulk_enrich_bar = ggpar(bulk_enrich_bar, legend = "right")


#put volcano and pathway enrich next to each other
figure = grid.arrange(bulk_volcano, bulk_enrich_bar, ncol = 2, widths = c(1, 1))
annotate_figure(figure, top = text_grob("Summary of DGE results in aggressive vs slowly proliferating primary tumours", size = 20, face = "bold"))



