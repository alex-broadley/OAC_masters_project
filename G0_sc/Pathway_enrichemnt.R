setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

library(escape)
library(pheatmap)
library(stringr)
library(ggplot2)
library(ggpubr)

#tum_categories = read.csv('G0_scored_tumorCells.csv')
tum_categories = read.csv('G0_scored_stringent.csv',row.names=1)
epi_counts = read.csv("epi_raw_counts_matrix.csv")

#make row counts gene names for both dataframes
rownames(tum_categories) = tum_categories$Sample
rownames(epi_counts) = epi_counts$X
epi_counts <- subset(epi_counts, select = -X)

#subset epithelial to get just tumor cell counts
tumor_counts = epi_counts[,which(colnames(epi_counts) %in% tum_categories$Sample)]
#write.csv(tumor_counts,'tum_counts.csv', row.names=TRUE)

#get matrix of the fast and G0 arrested cells
fast_G0_counts = tumor_counts[,which(colnames(tumor_counts) %in% tum_categories$Sample[which(tum_categories$strict_cycle_status %in% c('fast cycling','G0 arrested'))])]
fast_G0_counts = as.matrix(fast_G0_counts)
#write.csv(fast_G0_counts,'fast_G0_counts.csv', row.names=TRUE)


# ----------- MigSig Hallmark gene sets ---------------

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
G0_fast_metadata$Cycle_type = gsub("G0 arrested", "G0_arrested", G0_fast_metadata$Cycle_type)
G0_fast_metadata$Cycle_type = gsub("fast cycling", "fast_cycling", G0_fast_metadata$Cycle_type)
anno <- as.data.frame(G0_fast_metadata)
anno <- subset(anno, select = -Sample)
# define annotation colors

annoCol<-list(Cycle_type=c(G0_arrested='#286837', fast_cycling='#772989'))


pheatmap(t(enrichment.scores), 
         cluster_rows = FALSE,
         cluster_columns = FALSE,
         fontsize = 7,
         annotation_col = anno,
         show_colnames = FALSE,
         annotation_colors = annoCol,
         color = colorRampPalette(c("#16245a","white","#BD1E1E"))(100)) 


# ----------- Gene ontology (GO) terms -----------

all_GO_sets <- getGeneSets(library="C5", subcategory = "BP")

GO_selected_sets = c("GOBP-POSITIVE-REGULATION-OF-CELL-CYCLE-G1-S-PHASE-TRANSITION", "GOBP-SIGNAL-TRANSDUCTION-BY-P53-CLASS-MEDIATOR",
                  "GOBP-POSITIVE-REGULATION-OF-CELL-CYCLE-G2-M-PHASE-TRANSITION", 
                  "GOBP-DNA-DAMAGE-RESPONSE-SIGNAL-TRANSDUCTION-BY-P53-CLASS-MEDIATOR-RESULTING-IN-CELL-CYCLE-ARREST")

GO_sets <- all_GO_sets[GO_selected_sets]


names(GO_sets) <- gsub('-', '_', names(GO_sets))
names(GO_sets) <- gsub('GOBP_', '', names(GO_sets))
names(GO_sets) <- c("PosReg_G1S", "SigTrans_P53", "PosReg_G2M", "P53_DNAdamage_G0arrest")

#calculate enrichment score 
GO_enrichment.scores <- escape.matrix(as.matrix(tumor_counts), 
                                   gene.sets = GO_sets, 
                                   min.size = 5)

#scale enrichment scores for better heatmap plot
GO_enrichment.scores <- GO_enrichment.scores - rowMeans(GO_enrichment.scores)

enrich_metadata = merge(GO_enrichment.scores, tum_categories, by = 0) 
names(enrich_metadata)[12] = "Cycle_type"

enrich_metadata$Cycle_type = gsub("G0 arrested", "G0_arrested", enrich_metadata$Cycle_type)
enrich_metadata$Cycle_type = gsub("fast cycling", "fast_cycling", enrich_metadata$Cycle_type)


enrich_metadata$Cycle_type = factor(enrich_metadata$Cycle_type , levels=c("fast_cycling", "cycling", "G0_arrested"))

my_comparisons = list(c("fast_cycling", "cycling"), c("cycling", 'G0_arrested'), c("fast_cycling", "G0_arrested"))

colnames(enrich_metadata)

compare_means(PosReg_G1S ~ Cycle_type,  data = enrich_metadata)
compare_means(TOP2A ~ Cycle_type,  data = tumor_counts)

ggplot(enrich_metadata, aes(x=Cycle_type, y=PosReg_G1S, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 2200) +
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_bw()

ggplot(enrich_metadata, aes(x=Cycle_type, y=SigTrans_P53, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 2200) +
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_bw()

ggplot(enrich_metadata, aes(x=Cycle_type, y=PosReg_G2M, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 2600) +
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_bw()

ggplot(enrich_metadata, aes(x=Cycle_type, y=P53_DNAdamage_G0arrest, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 2300) +
  scale_fill_manual(values=c("#b30000", "#9494b8", "#000080")) +
  theme_bw()





