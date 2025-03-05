setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

library(escape)
library(pheatmap)
library(stringr)

tum_categories = read.csv('G0_scored_tumorCells.csv')
epi_counts = read.csv("epi_raw_counts_matrix.csv")

#make row counts gene names for both dataframes
rownames(tum_categories) = tum_categories$Sample
rownames(epi_counts) = epi_counts$X
epi_counts <- subset(epi_counts, select = -X)

#subset epithelial to get just tumor cell counts
tumor_counts = epi_counts[,which(colnames(epi_counts) %in% tum_categories$Sample)]
write.csv(tumor_counts,'tum_counts.csv', row.names=TRUE)

#get matrix of the fast and G0 arrested cells
fast_G0_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$strict_cycle_status %in% c('fast','slow'))])]
fast_G0_counts = as.matrix(G0_counts)
write.csv(fast_G0_counts,'fast_G0_counts.csv', row.names=TRUE)

#built in default is MigSig database
hallmark_geneset <- getGeneSets(library = "H")

#subset same hallmarks as Cenk
selected_sets = c("HALLMARK-KRAS-SIGNALING-DN", "HALLMARK-MYC-TARGETS-V2", 'HALLMARK-OXIDATIVE-PHOSPHORYLATION', "HALLMARK-DNA-REPAIR", 
                  "HALLMARK-UNFOLDED-PROTEIN-RESPONSE", "HALLMARK-MYC-TARGETS-V1", "ALLMARK-MTORC1-SIGNALING", "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY",
                  "HALLMARK-PI3K-AKT-MTOR-SIGNALING", "HALLMARK-PROTEIN-SECRETION", "HALLMARK-P53-PATHWAY", "HALLMARK-MITOTIC-SPINDLE", "HALLMARK-E2F-TARGETS", "HALLMARK-G2M-CHECKPOINT")

hallmark_geneset <- hallmark_geneset[selected_sets]

#remove NA gene set
hallmark_geneset <- hallmark_geneset[-7]


#format hallmark gene set names for nicer plotting
names(hallmark_geneset) <- gsub('HALLMARK', '', names(hallmark_geneset))
names(hallmark_geneset) <- gsub('-', ' ', names(hallmark_geneset))
names(hallmark_geneset) <- sub('.', '', names(hallmark_geneset))
names(hallmark_geneset)[8] <- "PI3K/AKT/MTOR SIGNALING"
names(hallmark_geneset) <- tolower(names(hallmark_geneset))
names(hallmark_geneset) <- str_to_title(names(hallmark_geneset))

#calculate enrichment score 
enrichment.scores <- escape.matrix(fast_G0_counts, 
                                   gene.sets = hallmark_geneset, 
                                   min.size = 5)

#scale enrichment scores for better heatmap plot
enrichment.scores <- enrichment.scores - rowMeans(enrichment.scores)

#create annotation dataframe to label heatmap
G0_fast_metadata = tum_categories[which(tum_categories$Sample %in% colnames(fast_G0_counts)), c('Sample','strict_cycle_status')]
names(G0_fast_metadata)[2] = 'Cycle_type'
anno <- as.data.frame(G0_fast_metadata)
anno <- subset(anno, select = -Sample)
# define annotation colors

annoCol<-list(Cycle_type=c(fast='#286837', slow='#772989'))


pheatmap(t(enrichment.scores), 
         cluster_rows = FALSE,
         cluster_columns = FALSE,
         fontsize = 7,
         annotation_col = anno,
         show_colnames = FALSE,
         annotation_colors = annoCol,
         color = colorRampPalette(c("#16245a","white","#BD1E1E"))(100)) 
