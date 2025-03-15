#sc DECAF on cancer cells by proliferation rate

library(scDECAF)
library(Seurat)
library(SingleCellExperiment)
library(escape)
library(scater)

setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

#G0 metadata for tumor cells
tum_categories = read.csv('G0_scored_tumorCells.csv',row.names=1)
#rows need to be genes, columns cells
tumor_counts = read.csv('tum_counts.csv', row.names = 1)


print(tum_categories[!(rownames(tum_categories) %in% colnames(tumor_counts)),])

#temp fix before rectify file production
tum_categories = tum_categories[(rownames(tum_categories) %in% colnames(tumor_counts)),]

#normalize counts -> should probably take these values from overall scanpy anndata object
tumor_counts = NormalizeData(tumor_counts)  

#create SCE object with metadata
tum_sce <- SingleCellExperiment(list(logcounts=as(as.matrix(tumor_counts), "dgCMatrix")), colData=DataFrame(tum_categories))

#load MigSig hallmark genesets - subset same hallmarks as Cenk
hallmark_geneset <- getGeneSets(library = "H")
selected_sets = c("HALLMARK-P53-PATHWAY", "HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION")
hallmark_geneset <- hallmark_geneset[selected_sets]

#take just the genes measured that are also in the genesets
filtered_tum_sce <- tum_sce[intersect(rownames(tum_sce), unique(as.vector(unlist(hallmark_geneset)))),]

#get HVGs - Don't know if should do this as just want all the genes from the geneset?
tum_varModel <-  scran::modelGeneVar(filtered_tum_sce, block=filtered_tum_sce$strict_cycle_status, density.weights=FALSE)
features_lst <- scran::getTopHVGs(tum_varModel , n=3000)

#cell embeddings - UMAP
filtered_tum_sce <- scater::runUMAP(filtered_tum_sce, scale=FALSE,  exprs_values = "logcounts")
cell_embed = as.data.frame(reducedDim(filtered_tum_sce, "UMAP"))




### --- Run scDECAF
data_df <- as.matrix(logcounts(filtered_tum_sce))
target <- scDECAF::genesets2ids(data_df[match(features_lst, rownames(data_df)),], hallmark_geneset) 
genesetThresh=5
standardise=TRUE
genesetDrop <- unlist(lapply(as.vector(colSums(target)), function(x){ifelse(x<genesetThresh, FALSE, TRUE)}))
target <- target[,genesetDrop]

output <- scDECAF::scDECAF(data = data_df, gs = target, standardize = standardise, 
                           hvg = features_lst, k = 10, embedding = cell_embed, cca.k=10,
                           n_components = 10, max_iter = 2, thresh = 0.5)

tum_sce_DECAF_results <- attr(output,"raw_scores")


selected_gs <- scDECAF::pruneGenesets(data = tum_sce, genesetlist = hallmark_geneset, hvg = rownames(data),
                             embedding = cell_embedding, min_gs_size = 3, lambda = exp(-3))




