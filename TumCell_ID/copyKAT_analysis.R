library(copykat)



setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

# copyKAT requires raw count expression matrix: rownames = GENE ids, colnames = cell ids
raw_counts = read.csv("epi_und_raw_counts.csv", row.names='X')
raw_counts_TME = read.csv("epi_und_TME_raw_counts.csv", row.names = 'X')

#loading SCEVAN results

SCEVAN_results_noTME = read.csv("SCEVAN_results_epi_und_noRef.csv", row.names = 'X')
SCEVAN_results_TME = read.csv("SCEVAN_results_epi_und_TME.csv", row.names = 'X')


### Running copyKAT -----------------------------------------------------------
# here run copyKAT with default parameters recommended on GitHUB

copykat.results_noTME <- copykat(rawmat=raw_counts, id.type="S", ngene.chr=5, 
                        win.size=25, KS.cut=0.1, sam.name="test", 
                        distance="euclidean", norm.cell.names="",output.seg="FLASE", 
                        plot.genes="TRUE", genome="hg20",n.cores=1)

copykat.results_TME <- copykat(rawmat=raw_counts_TME, id.type="S", ngene.chr=5, 
                                 win.size=25, KS.cut=0.1, sam.name="test", 
                                 distance="euclidean", norm.cell.names="",output.seg="FLASE", 
                                 plot.genes="TRUE", genome="hg20",n.cores=1)

### Analysing results -----------------------------------------------------------

# Benchmarking study preferred SCEVAN for break-point analysis, so most interested in tumor cell classification
# Aneuploid = tumor, diploid = normal

predictions_noTME <- data.frame(copykat.results_noTME$prediction)
predictions_noTME <- predictions_noTME[which(predictions_noTME$copykat.pred %in% c("aneuploid","diploid")),]  
copyKAT_tumorIDs_noTME <- predictions_noTME$cell.names[which(predictions_noTME$copykat.pred=="aneuploid")]

predictions_TME <- data.frame(copykat.results_TME$prediction)
predictions_TME <- predictions_TME[which(predictions_TME$copykat.pred %in% c("aneuploid","diploid")),]  
copyKAT_tumorIDs_TME <- predictions_TME$cell.names[which(predictions_TME$copykat.pred=="aneuploid")]


### Export results -----------------------------------------------------------

write.csv(copyKAT_tumorIDs_noTME, 'copyKAT_tumors_noTME.csv')
write.csv(copyKAT_tumorIDs_TME, 'copyKAT_tumors_TME.csv')



### Comparing results of SCEVAN and copyKAT in classification -----------------------------------------------------------


SCEVAN_tumorIDs_noTME = rownames(SCEVAN_results_noTME[which(SCEVAN_results_noTME$class == 'tumor'), ])
SCEVAN_tumorIDs_TME = rownames(SCEVAN_results_TME[which(SCEVAN_results_TME$class == 'tumor'), ])

#compare SCEVAN and copyKAT results w/o TME cells included in dataset
table(copyKAT_tumorIDs_noTME %in% SCEVAN_tumorIDs_noTME)
table(SCEVAN_tumorIDs_noTME %in% copyKAT_tumorIDs_noTME)

#compare SCEVAN and copyKAT results with TME cells included in dataset
table(copyKAT_tumorIDs_TME %in% SCEVAN_tumorIDs_TME)
table(SCEVAN_tumorIDs_TME %in% copyKAT_tumorIDs_TME)

#compare copyKAT results w/ and w/o TME
table(copyKAT_tumorIDs_TME %in% copyKAT_tumorIDs_noTME)
table(copyKAT_tumorIDs_noTME %in% copyKAT_tumorIDs_TME)

#compare SCEVAN results w/ and w/o TME
table(SCEVAN_tumorIDs_noTME %in% SCEVAN_tumorIDs_TME)
table(SCEVAN_tumorIDs_TME %in% SCEVAN_tumorIDs_noTME)




