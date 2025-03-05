library("dndscv")

#load datasets

load("barretts.vcf.RData")
load("mets.vcf.RData")
load("primaries.vcf.RData")

pooled_data = rbind(barretts_clean, primaries_clean, mets_clean)

####Getting VCF file names to reconstruct originals for SSB selection tools
barretts.vcf$Sample <- paste0(barretts.vcf$Sample, ".vcf")
mets.vcf$Sample <- paste0(mets.vcf$Sample, ".vcf")
primaries.vcf$Sample <- paste0(primaries.vcf$Sample, ".vcf")

write.csv(unique(barretts.vcf$Sample), 'B_sampleNames.csv') 
write.csv(unique(mets.vcf$Sample), 'M_sampleNames.csv') 
write.csv(unique(primaries.vcf$Sample), 'P_sampleNames.csv') 

####Getting stats on number of samples
barretts.vcf$split_ID = barretts.vcf$Sample
primaries.vcf$split_ID = primaries.vcf$Sample
mets.vcf$split_ID = mets.vcf$Sample

barretts.vcf$split_ID <- sapply(barretts.vcf$split_ID, function(x) strsplit(x,"_vs_")[[1]][1])
primaries.vcf$split_ID <- sapply(primaries.vcf$split_ID, function(x) strsplit(x,"_vs_")[[1]][1])
mets.vcf$split_ID <- sapply(mets.vcf$split_ID, function(x) strsplit(x,"_vs_")[[1]][1])

length(unique(barretts.vcf$Sample))
length(unique(primaries.vcf$Sample))
length(unique(mets.vcf$Sample))

##### Metastatic dataset

#extract only relevant rows
mets_clean = mets.vcf[c("#CHROM", "POS", "REF", "ALT", "Sample")]
#change column names
names(mets_clean) = c('chr', 'pos', 'ref', 'mut', 'sampleID')
#put the columns in the correct order as per the tutorial
mets_clean = mets_clean[c('sampleID', 'chr', 'pos','ref','mut')]
#correct chromosome name issue
mets_clean$chr <- gsub("chr","",as.vector(mets_clean$chr))

#run dndscv to find selected genes
dndsout_mets = dndscv(mets_clean, outp=3)

#extract gene names and associated statistics for all significantly positively selected genes
sel_cv_mets = dndsout_mets$sel_cv
signif_met_genes = sel_cv_mets[sel_cv_mets$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                "wmis_cv","wnon_cv","wspl_cv")]
sel_loc_mets = dndsout_mets$sel_loc
#remove indexing
rownames(signif_met_genes) = NULL

#save as RData file
save(signif_met_genes, file="signif_genes.metastases.RData")
save(sel_cv_mets, file = "sel_cv_mets.RData")
save(sel_loc_mets, file = "sel_loc_mets.RData")


##### Barretrs dataset


#extract only relevant rows
barretts_clean = barretts.vcf[c("#CHROM", "POS", "REF", "ALT", "Sample")]

#change column names
names(barretts_clean) = c('chr', 'pos', 'ref', 'mut', 'sampleID')

#put the columns in the correct order as per the tutorial
barretts_clean = barretts_clean[c('sampleID', 'chr', 'pos','ref','mut')]

#correct chromosome name issue
barretts_clean$chr <- gsub("chr","",as.vector(barretts_clean$chr))

#run dndscv to find selected genes
dndsout_barretts = dndscv(barretts_clean, outp=3)

#extract gene names and associated statistics for all significantly positively selected genes
sel_cv_barretts = dndsout_barretts$sel_cv
signif_barretts_genes = sel_cv_barretts[sel_cv_barretts$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                         "wmis_cv","wnon_cv","wspl_cv")]

sel_loc_barretts = dndsout_barretts$sel_loc
#remove indexing
rownames(signif_barretts_genes) = NULL

#save as RData file
save(signif_barretts_genes, file="signif_genes.barretts.RData")
save(sel_cv_barretts, file = "sel_cv_barretts.RData")
save(sel_loc_barretts, file = "sel_loc_barretts.RData")


##### Primary dataset

#extract only relevant rows and format issues
primaries_clean = primaries.vcf[c("#CHROM", "POS", "REF", "ALT", "Sample")]
#change column names
names(primaries_clean) = c('chr', 'pos', 'ref', 'mut', 'sampleID')
#put the columns in the correct order as per the tutorial
primaries_clean = primaries_clean[c('sampleID', 'chr', 'pos','ref','mut')]
#correct chromosome name issue
primaries_clean$chr <- gsub("chr","",as.vector(primaries_clean$chr))

#run dndscv to find selected genes
dndsout_primaries = dndscv(primaries_clean, outp=3)

#extract gene names and associated statistics for all significantly positively selected genes
sel_cv_primaries = dndsout_primaries$sel_cv
signif_primary_genes = sel_cv_primaries[sel_cv_primaries$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                  "wmis_cv","wnon_cv","wspl_cv")]

sel_loc_primaries = dndsout_primaries$sel_loc
#remove df indexing
rownames(signif_primary_genes) = NULL

#save as RData file
save(signif_primary_genes, file="signif_genes.primary.RData")
save(sel_cv_primaries, file ="sel_cv_primaries.RData")
save(sel_loc_primaries, file ="sel_loc_primaries.RData")


####extracting a table in which the genes significant in either barretts or primary 
#are included for both

#create a gene list with genes significant in either
extracted_gene_list = unique(c(sel_cv_primaries[sel_cv_primaries$qallsubs_cv<0.1, c("gene_name")], sel_cv_barretts[sel_cv_barretts$qallsubs_cv<0.1, c("gene_name")]))

#subset other dataframes and rename columns pre-merge
temp_primary = sel_cv_primaries[sel_cv_primaries$gene_name %in% extracted_gene_list, c("gene_name", "qallsubs_cv")]
names(temp_primary) = c('gene_name', 'primary_P')

temp_barretts = sel_cv_barretts[sel_cv_barretts$gene_name %in% extracted_gene_list, c("gene_name", "qallsubs_cv")]
names(temp_barretts) = c('gene_name', 'barretts_P')

#merge on gene name
signif_primary_barretts <- merge(temp_primary, temp_barretts,by="gene_name")

#save as Rdata object
save(signif_primary_barretts, file = 'signif_genes.primary_barretts.RData')


### Getting background genes for DAVID enrichment analysis
write.csv(sel_cv_mets$gene_name, 'all_metastatic_genes.csv') 

#getting foreground genes
write.csv(signif_barretts_genes$gene_name, 'sig_BE_genes.csv') 
write.csv(signif_primary_genes$gene_name, 'sig_prim_genes.csv') 
write.csv(signif_met_genes$gene_name, 'sig_met_genes.csv') 
write.csv(signif_pooled_genes$gene_name, 'sig_pooled_genes.csv') 

###### Pooled analysis

pooled_data = rbind(barretts_clean, primaries_clean, mets_clean)

dndsout_pooled = dndscv(pooled_data)
sel_cv_pooled = dndsout_pooled$sel_cv
signif_pooled_genes = sel_cv_pooled[sel_cv_pooled$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                           "wmis_cv","wnon_cv","wspl_cv")]
save(signif_pooled_genes, file = 'signif_pooled_genes.RData')

