##### calculating G0 arrest scores for bulk tumours samples
#then splits cohorts based on proportion of arrested cells in tumour
library(readxl)

setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')

#Load expression matrices - mat.expr has expression for scoring, oac.expr used for matching ID label
load("oac.expr.final.RData")
load("mat.expr.RData")

#write in non-Rdata format
write.csv(mat.expr, "/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data/BULK_TUM_RNA.csv")

#remove gene names from expression matrix, set as rownames instead
mat.expr.oac <- mat.expr[,-1]
mat.expr.oac <- as.matrix(mat.expr.oac)
rownames(mat.expr.oac) <- mat.expr$DNAid

# Read the G0 arrest signature:
G0sig <- read_xlsx("G0arrest_signature.xlsx")
G0sig.up <- G0sig[which(G0sig$Expression=="upregulated"),]$Gene
G0sig.down <- G0sig[which(G0sig$Expression=="downregulated"),]$Gene

# List of signatures:
sigs <- list()
sigs[[1]] <- G0sig.up
sigs[[2]] <- G0sig.down


# Calculate scores:
library(GSVA)
## build GSVA parameter object
gsvapar <- gsvaParam(t(mat.expr.oac), sigs, maxDiff=TRUE)



#check gene signatures match those in expr matrix
table(rownames(t(mat.expr.oac)) %in% sigs[[1]])
table(rownames(t(mat.expr.oac)) %in% sigs[[2]])

## calculate scores
gsva_es <- gsva(gsvapar)
rownames(gsva_es) <- c("G0.up","G0.down")
mat.scores <- t(gsva_es)
mat.scores <- data.frame(mat.scores)
mat.scores$Sample <- rownames(mat.scores)
mat.scores$ProliferativeCapacity <- mat.scores$G0.down-mat.scores$G0.up
# Prolif capacity - the higher the more proliferative, the lower the slower dividing;
# above 0 should be relatively fast/aggressive tumour
mat.scores$prolif_z = scale(mat.scores$ProliferativeCapacity, center = T, scale = T)

##### Splitting primary tumour cohorts

#load in vcf data
load("primaries.vcf.RData")
load("mets.vcf.RData")
load("barretts.vcf.RData")

### Map sample ID's to DNAid's used in VCF files
mapping <- unique(oac.expr[,c("Sample","DNAid")])
mat.scores.mapped <- merge(mat.scores, mapping, by.x="Sample", by.y="Sample", all.x=F, all.y=F)

# Split cohort by proliferative score - checked if any zero values and found none
slow_grow = mat.scores.mapped[mat.scores.mapped$ProliferativeCapacity < 0,]
quick_grow = mat.scores.mapped[mat.scores.mapped$ProliferativeCapacity > 0,]
zero_grow = mat.scores.mapped[mat.scores.mapped$ProliferativeCapacity == 0,]

# #Split the sample names in the vcf file by '_vs_' to allow matching

#making a column for the cancer sample ID and healthy reference tissue ID
primaries.vcf$split_ID = primaries.vcf$Sample
primaries.vcf$split_ID2 = primaries.vcf$Sample
primaries.vcf$split_ID <- sapply(primaries.vcf$split_ID, function(x) strsplit(x,"_vs_")[[1]][1])
primaries.vcf$split_ID2 <- sapply(primaries.vcf$split_ID, function(x) strsplit(x,"_vs_")[[1]][2])

### find the vcf rows with matching sample IDs
quick_grow_samples = primaries.vcf[primaries.vcf$split_ID %in% quick_grow$DNAid,]
slow_grow_samples = primaries.vcf[primaries.vcf$split_ID %in% slow_grow$DNAid,]

# Checking that no ID's match up with the healthy reference tissue instead
test_match = primaries.vcf[primaries.vcf$split_ID2 %in% quick_grow$DNAid,]
test_match = primaries.vcf[primaries.vcf$split_ID2 %in% slow_grow$DNAid,]

# checking there aren't any non-primary tumor matches 
test_match = barretts.vcf[barretts.vcf$split_ID %in% quick_grow$DNAid,]
test_match = barretts.vcf[barretts.vcf$split_ID %in% slow_grow$DNAid,]

test_match = mets.vcf[mets.vcf$split_ID %in% quick_grow$DNAid,]
test_match = mets.vcf[mets.vcf$split_ID %in% slow_grow$DNAid,]

#remove healthy ID column (no longer need)
quick_grow_samples = quick_grow_samples[ , !(names(quick_grow_samples) %in% c('split_ID2'))]
slow_grow_samples = slow_grow_samples[ , !(names(slow_grow_samples) %in% c('split_ID2'))]

#save split vcfs for use in selection analysis
save(quick_grow_samples, file = 'quick_grow_primaries.vcf.RData')
save(slow_grow_samples, file = 'slow_grow_primaries.vcf.RData')


#GSVA not working again -> get slow/fast growing sample IDs from these files
load('quick_grow_primaries.vcf.RData')
load('slow_grow_primaries.vcf.RData')

#classify tumors as fast/slow growing by splitting on score of 0 
#- here used whether or not in split data file as GSVA broke
expr_df = as.data.frame(mat.expr)
expr_df$proliferation <- ifelse(expr_df$DNAid %in% slow_grow_samples$split_ID, "slow", 
                           ifelse(expr_df$DNAid %in% quick_grow_samples$split_ID, "fast", 
                                  NA))


setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

#write classified tumors for downstream analysis
write.csv(expr_df, 'bulk_RNA_classified.csv', row.names = TRUE)


####################### GSVA broke - for figure scored in python and imported

bulk_scores = read.csv('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data/BULK_PROLIF_SCORE.csv', row.names = 'X')

names(bulk_scores)[1] = 'prolif_z'
bulk_scores$prolif_z_scaled = scale(bulk_scores$prolif_z, center = T, scale = T)

bulk_scores$Cycle_type = ""

#classify based on score of 0
for (i in 1:nrow(bulk_scores)){
  if (bulk_scores$prolif_z[i] >= 0) {
    bulk_scores$Cycle_type[i] = 'Aggressive'
  } else if (bulk_scores$prolif_z[i] <= 0) {
    bulk_scores$Cycle_type[i] = 'Slowly proliferating'
  }
}


# Plot histogram using same color scheme and theme as other figures
ggplot(data = bulk_scores, aes(prolif_z_scaled, fill = Cycle_type)) +
  scale_color_manual(values=c("#e60039", "#3366ff"))+
  scale_fill_manual(values=c("#e60039", "#3366ff")) +
  geom_histogram(breaks = c(seq(-2,3.5, 0.25)),color="black") +
  scale_x_continuous(breaks = -2:3.5) +
  theme_classic() +
  theme(axis.text=element_text(size=15), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  labs(
    x = "Proliferative score",
    y = "Counts",
    title = "Primary tumor proliferative score distribution",
    fill = 'Proliferation type'
  )

