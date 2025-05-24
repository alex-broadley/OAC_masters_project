#### Plots for significantly selected genes

setwd("/Users/alex/Documents/University/Year 4/BIOL0041/Code/BIOL0041")

#loading in RData findings
load("signif_genes.barretts.RData")
load("signif_genes.primary.RData")
load("signif_genes.metastases.RData")
load("signif_pooled_genes.RData")
signif_quick_genes = read.csv("signif_quick_genes.csv", row.names = 'X')
signif_slow_genes = read.csv("signif_slow_genes.csv", row.names = 'X')[,2:10]

##initialize libraries
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(ggpubr)

#take -log10 of p values for plotting
signif_barretts_genes$minusLogQall = -log10(signif_barretts_genes$qallsubs_cv)
signif_met_genes$minusLogQall = -log10(signif_met_genes$qallsubs_cv)
signif_primary_genes$minusLogQall = -log10(signif_primary_genes$qallsubs_cv)
signif_pooled_genes$minusLogQall = -log10(signif_pooled_genes$qallsubs_cv)
signif_quick_genes$minusLogQall = -log10(signif_quick_genes$qallsubs_cv)
signif_slow_genes$minusLogQall = -log10(signif_slow_genes$qallsubs_cv)

#fix issue where pvalues of zero where infinite to make plotting more manageable
signif_barretts_genes$minusLogQall = as.numeric(gsub(Inf, 9, as.vector(signif_barretts_genes$minusLogQall)))
signif_met_genes$minusLogQall = as.numeric(gsub(Inf, 13, as.vector(signif_met_genes$minusLogQall)))
signif_primary_genes$minusLogQall = as.numeric(gsub(Inf, 9, as.vector(signif_primary_genes$minusLogQall)))
signif_quick_genes$minusLogQall = as.numeric(gsub(Inf, 5, as.vector(signif_quick_genes$minusLogQall)))
signif_slow_genes$minusLogQall = as.numeric(gsub(Inf, 12, as.vector(signif_slow_genes$minusLogQall)))

#make sure values in these columns are numeric -> had a random seemingly one off issue where they weren't
signif_barretts_genes$minusLogQall = as.numeric(signif_barretts_genes$minusLogQall)
signif_met_genes$minusLogQall = as.numeric(signif_met_genes$minusLogQall)
signif_primary_genes$minusLogQall = as.numeric(signif_primary_genes$minusLogQall)
signif_pooled_genes$minusLogQall = as.numeric(signif_pooled_genes$minusLogQall)
signif_quick_genes$minusLogQall = as.numeric(signif_quick_genes$minusLogQall)
signif_slow_genes$minusLogQall = as.numeric(signif_slow_genes$minusLogQall)



pdf('scatter-signif-barrets.pdf')
dev.off()

ggscatter(signif_barretts_genes, 
          x = "wmis_cv", y = "minusLogQall", 
          xlab = "Strength of missense positive selection", ylab = "-log10(P-value)",
          label = "gene_name", repel = TRUE) + scale_y_continuous(breaks = seq(0, 11, by = 1))

ggscatter(signif_primary_genes, 
          x = "wmis_cv", y = "minusLogQall", 
          label = "gene_name", repel = TRUE) + scale_y_continuous(breaks = seq(0, 11, by = 1))

ggscatter(signif_met_genes[0:20,], 
          x = "wmis_cv", y = "minusLogQall", 
          label = "gene_name", repel = TRUE, size = 1) + scale_y_continuous(breaks = seq(0, 11, by = 1))

##### Plot ideas

# plotting points by their highest selection strength mutation type and coloring accordingly

most_selected_type_plot <- function(dfname, point_size, font_size, Title) {
  #initialise most_selected_score and type columns in the df so can edit in for loop
  dfname$most_selected_score = 0
  dfname$most_selected_type = 0
  
  #rename columns for nicer legend
  names(dfname)[7:9] = c('miss', 'non', 'splice')
  
  for (i in 1:nrow(dfname)) {
    #initialise composite label variable- used for creating legend/labelling points
    composite_label = ""
    #find which type of mutation is most strongly selected
    dfname$most_selected_score[i] = max(c(dfname$missense[i], dfname$nonsense[i], dfname$splice[i]))
    #take subset of the dataframe for the calculations - fixes indexing issue when using 'which'
    df_for_calcs = dfname[,c('miss', 'non', 'splice')]
    #create label showing all the mutation types within (10 for now) of the highest selection score
    if (length(names(df_for_calcs[,c('miss', 'non', 'splice')])[which(df_for_calcs[i, ] > 1)]) >= 1){
      highest_selection_mutations = names(df_for_calcs[,c('miss', 'non', 'splice')])[which(df_for_calcs[i, ] > 1)]
    }
    else {
      highest_selection_mutations = "None"
    }
    composite_label = paste(highest_selection_mutations, collapse = "-")
    #adding this label to the row
    dfname$most_selected_type[i] = composite_label
  }
  
  #make sure the label is a factor to allow plotting
  dfname$most_selected_type = as.factor(dfname$most_selected_type)
  
  print(dfname$most_selected_type)
  
  #create scatter plot, horizontal line indicates the p_value at which we decided the selection was significant
  plot = ggscatter(dfname, 
            x = "most_selected_score", y = "minusLogQall", 
            xlab = "Highest MLE dN/dS ratio across mutation types", ylab = "-log10(P-value)",
            main = Title,
            label = "gene_name", repel = TRUE, 
            color = "most_selected_type", 
            size = point_size, 
            font.label = c(font_size, "plain"),
            palette = c('red', 'blue', 'black')) + 
            scale_y_continuous(breaks = seq(0, 11, by = 1)) +
            scale_color_discrete(name = 'Mut. types with dN/dS > 1') +
            geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.3) +
            theme(axis.text=element_text(size=15), 
                  axis.text.x = element_text(vjust = 0.9, hjust=1),
                  axis.text.y = element_text(size=15),
                  axis.title.y = element_text(size=15),
                  axis.title.x = element_text(size = 15),
                  plot.title = element_text(hjust = 0.5, size = 15),
                  legend.text = element_text(size = 14))
  
  ggpar(plot, legend = "top")
} 


most_selected_type_plot <- function(dfname, point_size, font_size, Title) {
  #initialise most_selected_score and type columns in the df so can edit in for loop
  dfname$most_selected_score = 0
  dfname$most_selected_type = 0
  
  #rename columns for nicer legend
  names(dfname)[7:9] = c('miss', 'non', 'splice')
  
  for (i in 1:nrow(dfname)) {
    #initialise composite label variable- used for creating legend/labelling points
    composite_label = ""
    #find which type of mutation is most strongly selected
    dfname$most_selected_score[i] = max(c(dfname$missense[i], dfname$nonsense[i], dfname$splice[i]))
    #take subset of the dataframe for the calculations - fixes indexing issue when using 'which'
    df_for_calcs = dfname[,c('miss', 'non', 'splice')]
    #create label showing all the mutation types within (10 for now) of the highest selection score
    if (length(names(df_for_calcs[,c('miss', 'non', 'splice')])[which(df_for_calcs[i, ] > 1)]) >= 1){
      highest_selection_mutations = names(df_for_calcs[,c('miss', 'non', 'splice')])[which(df_for_calcs[i, ] > 1)]
    }
    else {
      highest_selection_mutations = "None"
    }
    composite_label = paste(highest_selection_mutations, collapse = "-")
    #adding this label to the row
    dfname$most_selected_type[i] = composite_label
  }
  
  #make sure the label is a factor to allow plotting
  dfname$most_selected_type = as.factor(dfname$most_selected_type)
  
  print(dfname$most_selected_type)
  
  #create scatter plot, horizontal line indicates the p_value at which we decided the selection was significant
  plot = ggplot(dfname, aes(x = most_selected_score, y = minusLogQall, label = gene_name, color = most_selected_type)) +
    geom_text_repel(max.overlaps = 20, size = font_size, force = 4) +
    geom_point() +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.3) +
    theme_bw() +
    theme(axis.text=element_text(size=15), 
          axis.text.x = element_text(vjust = 0.9),
          axis.text.y = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.title.x = element_text(size = 15),
          plot.title = element_text(hjust = 0.5, size = 15),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)) + 
    labs(x = "Highest MLE dN/dS ratio across mutation types", 
         y = "-log10(P-value)", 
         color = 'Mut. types with dN/dS > 1',
         title = Title) +
    scale_colour_manual(values=c("#00BB38", "#619DFF", "darkgrey"))
  
  ggpar(plot, legend = "top")
} 

# Candidate oncogene/TSP plot ---- NOT USED FOR DISS

candidate_geneType_plot <- function(dfname, point_size, font_size, Title) {
  #intialise column of for candidate labels and most selected score
  dfname$gene_type = ""
  dfname$most_selected_score = 0
  print(dfname$gene_type)
  
  #rename columns for nicer legend
  names(dfname)[7:9] = c('missense', 'nonsense', 'splice')
  
  for (i in 1:nrow(dfname)) {
    #if positive selection on nonesense or splice mutations, label as candidate tumor suppressor
    if (dfname$nonsense[i] > 1 | dfname$splice[i] > 1){
      dfname$gene_type[i] = 'Tumor suppressor'
    }
    else if (dfname$missense[i] > 1) {
      dfname$gene_type[i] = 'Oncogene'
    }
    #get most selected mutation type score for x axis
    dfname$most_selected_score[i] = max(c(dfname$missense[i], dfname$nonsense[i], dfname$splice[i]))
  }
  

  #make sure the label is a factor to allow plotting
  dfname$gene_type = as.factor(dfname$gene_type)
  
  print(dfname$gene_type)
  
  #create scatter plot, horizontal line indicates the p_value at which we decided the selection was significant
  plot = ggscatter(dfname, 
                   x = "most_selected_score", y = "minusLogQall", 
                   xlab = "Most strongly selected mutation type(s)", ylab = "-log10(P-value)",
                   main = Title,
                   label = "gene_name", repel = TRUE, color = "gene_type", size = point_size, font.label = c(font_size, "plain")) + 
    scale_y_continuous(breaks = seq(0, 11, by = 1)) +
    scale_color_discrete(name = 'Candidate gene type:') +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.3) +
    theme(plot.title = element_text(hjust = 0.5, size = 15))
  
  ggpar(plot, legend = "bottom")
} 


#running plots ---------------------------------------------------------------------

#candidate gene type - NOT USED IN DISS --------------------------------------

candidate_geneType_plot(signif_barretts_genes, 1, 12, "Barrett's Oesophagus")
candidate_geneType_plot(signif_primary_genes, 1, 12, "Prirmary OAC")
candidate_geneType_plot(signif_slow_genes, 1, 12, "Slowly proliferating tumors")
candidate_geneType_plot(signif_quick_genes, 1, 12, "Aggressively proliferating tumors")


#most selected type  ------------------------------------------------------------
mSelected_BO = most_selected_type_plot(signif_barretts_genes, 1, 5, "Mutation types under positive selection")
mSelected_pOAC = most_selected_type_plot(signif_primary_genes, 0.8, 3.8, "Mutation types under positive selection")
mSelected_met = most_selected_type_plot(signif_met_genes[0:30,], 0.8, 4.5,  "Mutation types under positive selection")
mSelected_quick = most_selected_type_plot(signif_quick_genes, 0.8, 5,  "Aggressively proliferating primary OAC")
mSelected_slow = most_selected_type_plot(signif_slow_genes, 0.8, 5,  "Slowly proliferating primary OAC")

#Plot to compare the genes significant in primary vs barretts -------------------
load("signif_primary_barretts.RData")

#convert to -log10
class(signif_primary_barretts$primary_P)
signif_primary_barretts$primaryLog = -log10(signif_primary_barretts$primary_P)
signif_primary_barretts$barrettsLog = -log10(signif_primary_barretts$barretts_P)

signif_primary_barretts$primaryLog = as.numeric(gsub(Inf, 9, as.vector(signif_primary_barretts$primaryLog)))
signif_primary_barretts$barrettsLog = as.numeric(gsub(Inf, 9, as.vector(signif_primary_barretts$barrettsLog)))

pdf('singif_primary_barretts_compared.pdf')
ggscatter(signif_primary_barretts, 
          x = "barrettsLog", y = "primaryLog", 
          xlab = "Significance in barretts", ylab = "Significance in primary",
          label = 'gene_name', repel = TRUE, size = 0.6, font.label = c(8, "plain")) +
          geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.3) +
          geom_vline(xintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.3)
dev.off()

#### plots of number of mutations by gene in a given stage -----------------------------------------------------------------------------------------------------

mutation_number_plot = function(df, m_types, Title, xaxis_fs){
  df$m_number = 0
  for (i in 1:nrow(df)){
    df$m_number[i] = sum(df[i,m_types])
  }
  
  ggbarplot(df, x="gene_name",y='m_number',fill = 'minusLogQall', main = Title, xlab = "", ylab = "Total non-synonymous mutations in cohort", sort.val = 'desc') +
    scale_fill_gradient2(low = "steelblue", high = "darkred", name = "-log10(P-value)") +
    theme(axis.text=element_text(size=15), 
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size = xaxis_fs),
          axis.text.y = element_text(size=15),
          axis.title.y = element_text(size=15),
          plot.title = element_text(hjust = 0.5, size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))
}

BO_m_numbers = mutation_number_plot(signif_barretts_genes, m_types = c('n_mis', 'n_non', 'n_spl'), Title = "Mutation prevalence", xaxis_fs = 15)
prim_m_numbers = mutation_number_plot(signif_primary_genes, m_types = c('n_mis', 'n_non', 'n_spl'), Title = "Mutation prevalence", xaxis_fs = 10)
met_m_numbers = mutation_number_plot(signif_met_genes[0:30,], m_types = c('n_mis', 'n_non', 'n_spl'), Title = "Mutation prevalence", xaxis_fs = 13)

figure = grid.arrange(BO_m_numbers, prim_m_numbers, ncol = 2)

annotate_figure(figure,
                top = text_grob("Number of non-synonymous mutations in genes significantly selected in disease stage cohorts", size = 15, face = "bold"))

### Combining relevant figures into panels

## BO prevalence mutation type plot
figure = grid.arrange(BO_m_numbers, mSelected_BO, ncol = 2)
annotate_figure(figure, top = text_grob("Barrett's Oesophagus", size = 20, face = "bold"))

## Primary tumors prevalence mutation type plot
figure = grid.arrange(prim_m_numbers, mSelected_pOAC, ncol = 2)
annotate_figure(figure, top = text_grob("Primary OAC", size = 20, face = "bold"))

## Metastatic tumors prevalence mutation type plot
figure = grid.arrange(met_m_numbers, mSelected_met, ncol = 2)
annotate_figure(figure, top = text_grob("Metastatic OAC (30 most significantly selected genes)", size = 20, face = "bold"))

#Mutation type plots for slow and aggressively dividing tumors
figure = grid.arrange(mSelected_quick, mSelected_slow, ncol = 2)
annotate_figure(figure, top = text_grob("Genes and mutation types under positive selection in aggressively and slowly proliferating tumors", size = 20, face = "bold"))


