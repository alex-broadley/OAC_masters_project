#### Pathway enrichment analysis from output of pyDeseq2 DGE analysis 
library(pathfindR)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

DGE_results = read.csv("DGE_G0_fast_results.csv", row.names = 'X')

sig_DGE <- data.frame(DGE_results[which(DGE_results$padj<0.05),
                             c("log2FoldChange","padj")])

sig_DGE$Gene = rownames(sig_DGE)

sig_DGE$padj <- as.numeric(df$padj)
sig_DGE$log2FoldChange <- as.numeric(df$log2FoldChange)


#run pathdindR w/ KEGG pathways
output_df_KEGG <- run_pathfindR(sig_DGE[,c("Gene","log2FoldChange","padj")],
                           output_dir="infection_vs_control", gene_sets="KEGG")
od_KEGG <- data.frame(output_df)


#create term gene graph from the results
term_gene_graph(result_df = od, use_description = TRUE)
enrichment_chart(result_df = output_df, top_terms = 50)

visualize_terms()

#run pathdindR w/ Reactome pathways
output_df_Reactome <- run_pathfindR(sig_DGE[,c("Gene","log2FoldChange","padj")],
                                output_dir="infection_vs_control", gene_sets="Reactome")

od_Reactome <- data.frame(output_df)
