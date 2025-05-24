setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

library(escape)
library(pheatmap)
library(stringr)
library(ggplot2)
library(ggpubr)

#read in tumor cell proliferative labels and normalised expression scores
tum_categories = read.csv('G0_agreement_quartiles.csv',row.names=1)
tum_counts = read.csv("agreement_tum_SCRAN_counts.csv", row.names='X')

#check that counts are normalised
max(tum_counts)

#ensure only have tumor cells in expression dataframe
tum_counts <- tum_counts[,colnames(tum_counts) %in% rownames(tum_categories)]

#get matrix of the fast and G0 arrested cells only - (no 'Cycling')
fast_G0_counts = tum_counts[,which(colnames(tum_counts) %in% rownames(tum_categories)[which(tum_categories$Cycle_type %in% c('fast cycling','G0 arrested'))])]
fast_G0_counts = as.matrix(fast_G0_counts)

# ----------- MigSig Hallmark gene sets ---------------

#get genes sets - built in default is MigSig database
hallmark_geneset <- getGeneSets(library = "H")

#subset same hallmarks as Cenk breat cancer dormancy paper
selected_sets = c("HALLMARK-KRAS-SIGNALING-DN", "HALLMARK-MYC-TARGETS-V2", 'HALLMARK-OXIDATIVE-PHOSPHORYLATION', "HALLMARK-DNA-REPAIR", 
                  "HALLMARK-UNFOLDED-PROTEIN-RESPONSE", "HALLMARK-MYC-TARGETS-V1", "ALLMARK-MTORC1-SIGNALING", "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY",
                  "HALLMARK-PI3K-AKT-MTOR-SIGNALING", "HALLMARK-PROTEIN-SECRETION", "HALLMARK-P53-PATHWAY", "HALLMARK-MITOTIC-SPINDLE", "HALLMARK-E2F-TARGETS", "HALLMARK-G2M-CHECKPOINT")

#subset gene sets to those selected
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

#calculate enrichment score  w/ ssGSEA
enrichment.scores <- escape.matrix(fast_G0_counts, 
                                   method = 'ssGSEA',
                                   gene.sets = hallmark_geneset,
                                   normalize = FALSE,
                                   min.size = 5)

#scale enrichment scores for better heatmap plot
enrichment.scores <- scale(enrichment.scores)

#create annotation dataframe to label heatmap
tum_categories$Sample = rownames(tum_categories)

G0_fast_metadata = tum_categories[which(rownames(tum_categories) %in% colnames(fast_G0_counts)), c('Sample','Cycle_type')]
G0_fast_metadata$Cycle_type = gsub("G0 arrested", "G0_arrested", G0_fast_metadata$Cycle_type)
G0_fast_metadata$Cycle_type = gsub("fast cycling", "fast_cycling", G0_fast_metadata$Cycle_type)
anno <- as.data.frame(G0_fast_metadata)
anno <- subset(anno, select = -Sample)
# define annotation colors

annoCol<-list(Cycle_type=c(G0_arrested='#3366ff', fast_cycling='#e60039'))

#plot heatmap of these scores - scale rows as recommended for ssGSEA
pheatmap(t(enrichment.scores), 
         cluster_rows = FALSE,
         cluster_columns = TRUE,
         scale = 'row',
         fontsize = 10,
         annotation_col = anno,
         annotation_legend = TRUE,
         show_colnames = FALSE,
         annotation_colors = annoCol,
         color = colorRampPalette(c("#16245a","white","#BD1E1E"))(100),
         legend_labels = c('G0 arrested', 'Fast cycling')) 


# ----------- Gene ontology (GO) terms -----------

#get gene ontology gene sets for biological processes
all_GO_sets <- getGeneSets(library="C5", subcategory = "BP")

#select same gene sets as Cenk paper again
GO_selected_sets = c("GOBP-POSITIVE-REGULATION-OF-CELL-CYCLE-G1-S-PHASE-TRANSITION", "GOBP-SIGNAL-TRANSDUCTION-BY-P53-CLASS-MEDIATOR",
                  "GOBP-POSITIVE-REGULATION-OF-CELL-CYCLE-G2-M-PHASE-TRANSITION", 
                  "GOBP-DNA-DAMAGE-RESPONSE-SIGNAL-TRANSDUCTION-BY-P53-CLASS-MEDIATOR-RESULTING-IN-CELL-CYCLE-ARREST")
#subset to selected sets
GO_sets <- all_GO_sets[GO_selected_sets]

#rename for nicer plotting
names(GO_sets) <- gsub('-', '_', names(GO_sets))
names(GO_sets) <- gsub('GOBP_', '', names(GO_sets))
names(GO_sets) <- c("PosReg_G1S", "SigTrans_P53", "PosReg_G2M", "P53_DNAdamage_G0arrest")

#calculate enrichment score 
GO_enrichment.scores <- escape.matrix(as.matrix(tum_counts), 
                                   gene.sets = GO_sets, 
                                   method = 'ssGSEA',
                                   min.size = 5,
                                   normalize = FALSE)

#merge scores with proliferation labels
enrich_metadata = merge(GO_enrichment.scores, tum_categories, by = 0) 
names(enrich_metadata)[11] = "Cycle_type"

#rename as plotting issues without '_' in place of spaces
enrich_metadata$Cycle_type = gsub("G0 arrested", "G0_arrested", enrich_metadata$Cycle_type)
enrich_metadata$Cycle_type = gsub("fast cycling", "fast_cycling", enrich_metadata$Cycle_type)

#set factor levels for X axis order
enrich_metadata$Cycle_type = factor(enrich_metadata$Cycle_type , levels=c("fast_cycling", "cycling", "G0_arrested"))

#set up comparisons for kruskal-wallis overlays
my_comparisons = list(c("fast_cycling", "cycling"), c("cycling", 'G0_arrested'), c("fast_cycling", "G0_arrested"))

#create violin for G1/S gene set
posG1S = ggplot(enrich_metadata, aes(x=Cycle_type, y=PosReg_G1S, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8, alpha = 0.3) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 6.5) + 
  labs(x = "", y = "NES", title = "Positive regulation G1/S") + 
  scale_fill_manual(values=c("#e60039", "#c2c2d6", "#3366ff")) +
  theme_bw() +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  scale_x_discrete(breaks=c('fast_cycling', 'cycling', 'G0_arrested'),
                   labels=c("Fast cycling", "Cycling", "G0 arrested"))

#create violin for P53 sig transduction gene set
sigP53 = ggplot(enrich_metadata, aes(x=Cycle_type, y=SigTrans_P53, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8, alpha = 0.3) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 6.5) + 
  labs(x = "", y = "NES", title = "P53 signal transduction") + 
  scale_fill_manual(values=c("#e60039", "#c2c2d6", "#3366ff")) +
  theme_bw() + 
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  scale_x_discrete(breaks=c('fast_cycling', 'cycling', 'G0_arrested'),
                   labels=c("Fast cycling", "Cycling", "G0 arrested"))

#create violin for G2/M gene set
posG2M = ggplot(enrich_metadata, aes(x=Cycle_type, y=PosReg_G2M, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8, alpha = 0.3) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 6.5) + 
  labs(x = "", y = "NES", title = "Positive regulation G2/M") + 
  scale_fill_manual(values=c("#e60039", "#c2c2d6", "#3366ff")) +
  theme_bw() +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.position = "none") +
  scale_x_discrete(breaks=c('fast_cycling', 'cycling', 'G0_arrested'),
                   labels=c("Fast cycling", "Cycling", "G0 arrested"))

#create violin for P53 DNA damage arrest gene set
p53_damage = ggplot(enrich_metadata, aes(x=Cycle_type, y=P53_DNAdamage_G0arrest, fill = Cycle_type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="black", fill="white") +
  geom_point(position = position_jitter(seed = 1, width = 0.2), size = 0.8, alpha = 0.3) +
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 6.5) + 
  labs(x = "", y = "NES", title = "P53 DNA damage response") + 
  scale_fill_manual(values=c("#e60039", "#c2c2d6", "#3366ff")) +
  theme_bw() + 
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(family = "Helvetica", size = (15), hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.position = "none") +
    scale_x_discrete(breaks=c('fast_cycling', 'cycling', 'G0_arrested'),
                   labels=c("Fast cycling", "Cycling", "G0 arrested"))



library(grid)
# Put all four plots together in grid with two columns - have cell cycle on top row and p53 on bottom
figure = grid.arrange(posG1S, posG2M, sigP53, p53_damage, ncol = 2)





