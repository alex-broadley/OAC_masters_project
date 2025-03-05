## Retreives dataframe row for a gene of interest, used for exploring results 
# on command line
retrieve_gene_info = function(gene_name, df_name){
  print(df_name[df_name$gene_name %in% c(gene_name),])
}

#load in datasets - have not done
cosmic_genes = read.csv("Census_allThu Oct 24 12_28_26 2024.csv")
load("signif_genes.primary.RData")
load("signif_genes.metastases.RData")
load("signif_genes.barretts.RData")

# test to make sure that gene symbols are in suitable format for matching
control_not_cosmic = sel_cv_primaries[!(sel_cv_primaries$gene_name %in% cosmic_genes$Gene.Symbol), ]
not_removed_from_control = cosmic_genes[!(cosmic_genes$Gene.Symbol %in% sel_cv_primaries$gene_name), ]

write.csv(not_removed_from_control, 'cosmic_symbol_mismatch.csv')

#### Primary

#remove decimal points preventing some genes from being recognized
signif_primary_cosmic_removal = signif_primary_genes
signif_primary_cosmic_removal$gene_name <- gsub("\\..*", "", as.vector(signif_primary_cosmic_removal$gene_name))
#remove genes found in cosmic
primary_not_cosmic = signif_primary_cosmic_removal[!(signif_primary_cosmic_removal$gene_name %in% cosmic_genes$Gene.Symbol), ]
rownames(primary_not_cosmic) = NULL

write.csv(primary_not_cosmic, "primary_not_cosmic.csv")

##### Barretts

#remove decimal points preventing some genes from being recognized
signif_BE_cosmic_removal = signif_barretts_genes
signif_BE_cosmic_removal$gene_name <- gsub("\\..*", "", as.vector(signif_BE_cosmic_removal$gene_name))
#remove genes found in cosmic
BE_not_cosmic = signif_BE_cosmic_removal[!(signif_BE_cosmic_removal$gene_name %in% cosmic_genes$Gene.Symbol), ]
rownames(BE_not_cosmic) = NULL

write.csv(BE_not_cosmic, "BE_not_cosmic.csv")

##### Metastatic

#remove decimal points preventing some genes from being recognized
signif_met_cosmic_removal = signif_met_genes
signif_met_cosmic_removal$gene_name <- gsub("\\..*", "", as.vector(signif_met_cosmic_removal$gene_name))
#remove genes found in cosmic
met_not_cosmic = signif_met_cosmic_removal[!(signif_met_cosmic_removal$gene_name %in% cosmic_genes$Gene.Symbol), ]
rownames(met_not_cosmic) = NULL

write.csv(met_not_cosmic, "met_not_cosmic.csv")






