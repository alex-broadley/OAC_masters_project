##### Running dNdScv on bulk G0 stratified cohorts
library(dndscv)

setwd('/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041')

load('quick_grow_primaries.vcf.RData')
load('slow_grow_primaries.vcf.RData')

quick_clean = quick_grow_samples[c("#CHROM", "POS", "REF", "ALT", "Sample")]
slow_clean = slow_grow_samples[c("#CHROM", "POS", "REF", "ALT", "Sample")]

#correct column names
names(quick_clean) = c('chr', 'pos', 'ref', 'mut', 'sampleID')
names(slow_clean) = c('chr', 'pos', 'ref', 'mut', 'sampleID')

#putting columns in correct order
quick_clean = quick_clean[c('sampleID', 'chr', 'pos','ref','mut')]
slow_clean = slow_clean[c('sampleID', 'chr', 'pos','ref','mut')]

#reset column indexes
rownames(quick_clean) <- NULL
rownames(slow_clean) <- NULL

#correct chromosome labelling
quick_clean$chr <- gsub("chr","",as.vector(quick_clean$chr))
slow_clean$chr <- gsub("chr","",as.vector(slow_clean$chr))


### Run dndscv on the two cohorts ------------------------------------------------------
dndsout_quick = dndscv(quick_clean, outp=3)
dndsout_slow = dndscv(slow_clean, outp=3)

sel_cv_quick = dndsout_quick$sel_cv
signif_quick_genes = sel_cv_quick[sel_cv_quick$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                  "n_syn", "n_mis", "n_non", "n_spl",
                                                                  "wmis_cv","wnon_cv","wspl_cv")]

sel_cv_slow = dndsout_slow$sel_cv
signif_slow_genes = sel_cv_slow[sel_cv_slow$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                               "n_syn", "n_mis", "n_non", "n_spl",
                                                               "wmis_cv","wnon_cv","wspl_cv")]


### Output list of significantly selected genes
write.csv(signif_quick_genes, 'signif_quick_genes.csv') 
write.csv(signif_slow_genes, 'signif_slow_genes.csv')

