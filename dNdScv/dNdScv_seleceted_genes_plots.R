#### Plots for significantly selected genes

#loading in RData findings
load("signif_genes.barretts.RData")
load("signif_genes.primary.RData")
load("signif_genes.metastases.RData")
load("signif_pooled_genes.RData")

##initialize libraries
library(ggplot2)
library(ggpubr)

#### Calculate dN/dS for all mutation types together for each gene individually

calculate_total_dNdS = function(dndsout_table) {
  dndsout_table$total_dNdS = (dndsout_table$n_mis + dndsout_table$n_non + dndsout_table$n_spl) / dndsout_table$n_syn
  return(dndsout_table)
}

negative_selection_primaries = calculate_total_dNdS(dndsout_primaries$sel_cv)

negative_selection_primaries <- negative_selection_primaries[negative_selection_primaries$total_dNdS <0.1, ]
negative_selection_primaries <- negative_selection_primaries[negative_selection_primaries$total_dNdS >0, ]

#take -log10 of p values for plotting
signif_barretts_genes$minusLogQall = -log10(signif_barretts_genes$qallsubs_cv)
signif_met_genes$minusLogQall = -log10(signif_met_genes$qallsubs_cv)
signif_primary_genes$minusLogQall = -log10(signif_primary_genes$qallsubs_cv)
signif_pooled_genes$minusLogQall = -log10(signif_pooled_genes$qallsubs_cv)



#fix issue where pvalues of zero where infinite to make plotting more manageable
signif_barretts_genes$minusLogQall = as.numeric(gsub(Inf, 9, as.vector(signif_barretts_genes$minusLogQall)))
signif_met_genes$minusLogQall = as.numeric(gsub(Inf, 13, as.vector(signif_met_genes$minusLogQall)))
signif_primary_genes$minusLogQall = as.numeric(gsub(Inf, 9, as.vector(signif_primary_genes$minusLogQall)))
signif_pooled_genes$minusLogQall = as.numeric(gsub(Inf, 9, as.vector(signif_pooled_genes$minusLogQall)))

#make sure values in these columns are numeric -> had a random seemingly one off issue where they weren't
signif_barretts_genes$minusLogQall = as.numeric(signif_barretts_genes$minusLogQall)
signif_met_genes$minusLogQall = as.numeric(signif_met_genes$minusLogQall)
signif_primary_genes$minusLogQall = as.numeric(signif_primary_genes$minusLogQall)
signif_pooled_genes$minusLogQall = as.numeric(signif_pooled_genes$minusLogQall)



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

most_selected_type_plot <- function(dfname, point_size, font_size) {
  #initialise most_selected_score and type columns in the df so can edit in for loop
  dfname$most_selected_score = 0
  dfname$most_selected_type = 0
  
  #rename columns for nicer legend
  names(dfname)[3:5] = c('missense', 'nonsense', 'splice')
  
  for (i in 1:nrow(dfname)) {
    #initialise composite label variable- used for creating legend/labelling points
    composite_label = ""
    #find which type of mutation is most strongly selected
    dfname$most_selected_score[i] = max(c(dfname$missense[i], dfname$nonsense[i], dfname$splice[i]))
    #take subset of the dataframe for the calculations - fixes indexing issue when using 'which'
    df_for_calcs = dfname[,c('missense', 'nonsense', 'splice')]
    #create label showing all the mutation types within (10 for now) of the highest selection score
    highest_selection_mutations = names(df_for_calcs[,c('missense', 'nonsense', 'splice')])[which(df_for_calcs[i, ] >= (dfname$most_selected_score[i]-10))]
    composite_label = paste(highest_selection_mutations, collapse = "-")
    #assing this label to the row
    dfname$most_selected_type[i] = composite_label
  }
  
  #make sure the label is a factor to allow plotting
  dfname$most_selected_type = as.factor(dfname$most_selected_type)
  
  #create scatter plot, horizontal line indicates the p_value at which we decided the selection was significant
  ggscatter(dfname, 
            x = "most_selected_score", y = "minusLogQall", 
            xlab = "Strength of positive selection (for most highly selected mutation type)", ylab = "-log10(P-value)",
            label = "gene_name", repel = TRUE, color = "most_selected_type", size = point_size, font.label = c(font_size, "plain")) + 
            scale_y_continuous(breaks = seq(0, 11, by = 1)) +
            scale_color_discrete(name = 'Most highly selected mutation type(s)') +
            geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.3)
}

#running plots
most_selected_type_plot(signif_barretts_genes, 1, 10)
most_selected_type_plot(signif_primary_genes, 0.8, 10)
most_selected_type_plot(signif_met_genes[0:50,], 0.8, 10)
most_selected_type_plot(signif_pooled_genes, 0.8, 10)



#Plot to compare the genes significant in primary vs barretts and
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


#adding whether or not the genes significant in barrets or primaries are significant in metastatic samples to the plot
sel_cv_mets_subset = sel_cv_mets[, c('gene_name', 'qallsubs_cv')]
names(sel_cv_mets_subset) = c('gene_name', 'mets_P')
signif_primary_barretts <- merge(signif_primary_barretts, sel_cv_mets_subset,by="gene_name")

signif_primary_barretts$sig_met <- with(signif_primary_barretts, ifelse(mets_P < 0.1, 'Yes', 'No'))

ggscatter(signif_primary_barretts, 
          x = "barrettsLog", y = "primaryLog", 
          xlab = "Significance in BE", ylab = "Significance in primary tumour", color = 'sig_met',
          label = 'gene_name', repel = TRUE, size = 0.6, font.label = c(8, "plain")) +
          geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", linewidth = 0.3) + 
          geom_vline(xintercept = -log10(0.1), linetype = "dashed", color = "black", linewidth = 0.3) +
          scale_color_discrete(name = 'Significantly selected in metastatic samples', palette.colors(palette = 'okabe'))


#### plots of number of mutations by gene in a given stage

mutation_number_plot = function(df, m_types, Title){
  df$m_number = 0
  for (i in 1:nrow(df)){
    df$m_number[i] = sum(df[i,m_types])
  }
  
  ggbarplot(df, x="gene_name",y='m_number',fill = 'minusLogQall', main = Title, xlab = "Gene name", ylab = "Number of mutations in cohort", sort.val = 'desc') +
    scale_fill_gradient2(low = "steelblue", high = "darkred", name = "-log10(P-value)") +
    theme(axis.text=element_text(size=10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme(plot.title = element_text(hjust = 0.5))
    
  
}

mutation_number_plot(signif_barretts_genes, m_types = c('wnon_cv', 'wmis_cv', 'wspl_cv'), Title = "Mutation prevalence by gene in BE")
mutation_number_plot(signif_primary_genes, m_types = c('wnon_cv', 'wmis_cv', 'wspl_cv'), Title = "Mutation prevalence by gene in Primary Tumour")
mutation_number_plot(signif_met_genes[0:40,], m_types = c('wmis_cv', 'wnon_cv', 'wspl_cv'), Title = "Mutation prevalence by gene in Metastatic samples")

names(signif_primary_genes)

names(signif_barretts_genes)






