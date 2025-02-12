###########
### Try G0 arrest score:
library(readxl)

setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')

load("oac.expr.final.RData")
load("mat.expr.RData")

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

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar)
rownames(gsva_es) <- c("G0.up","G0.down")
mat.scores <- t(gsva_es)
mat.scores <- data.frame(mat.scores)
mat.scores$Sample <- rownames(mat.scores)
mat.scores$ProliferativeCapacity <- mat.scores$G0.down-mat.scores$G0.up
mat.scores$prolif_z = scale(mat.scores$ProliferativeCapacity, center = T, scale = T)
# Prolif capacity - the higher the more proliferative, the lower the slower dividing;
# above 0 should be relatively fast/aggressive tumour

hist(mat.scores$prolif_z)

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
nrow(slow_grow)
nrow(quick_grow)
nrow(zero_grow)

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

# checking there aren't any non-primary tumour matches 
test_match = barretts.vcf[barretts.vcf$split_ID %in% quick_grow$DNAid,]
test_match = barretts.vcf[barretts.vcf$split_ID %in% slow_grow$DNAid,]

test_match = mets.vcf[mets.vcf$split_ID %in% quick_grow$DNAid,]
test_match = mets.vcf[mets.vcf$split_ID %in% slow_grow$DNAid,]

#remove healthy ID column (no longer need)
quick_grow_samples = quick_grow_samples[ , !(names(quick_grow_samples) %in% c('split_ID2'))]
slow_grow_samples = slow_grow_samples[ , !(names(slow_grow_samples) %in% c('split_ID2'))]

save(quick_grow_samples, file = 'quick_grow_primaries.vcf.RData')
save(slow_grow_samples, file = 'slow_grow_primaries.vcf.RData')


