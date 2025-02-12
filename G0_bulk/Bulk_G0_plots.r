### BULK G0 arrest plots 
library(ggplot2)
library(ggpubr)

slow_genes = read.csv('signif_slow_genes.csv')
quick_genes = read.csv('signif_quick_genes.csv')

slow_genes$qallsubs_cv = -log10(slow_genes$qallsubs_cv)
slow_genes$qallsubs_cv = as.numeric(gsub(Inf, 9, as.vector(slow_genes$qallsubs_cv)))

quick_genes$qallsubs_cv = -log10(quick_genes$qallsubs_cv)
quick_genes$qallsubs_cv = as.numeric(gsub(Inf, 9, as.vector(quick_genes$qallsubs_cv)))

most_selected_type_plot <- function(dfname, point_size, font_size) {
  #initialise most_selected_score and type columns in the df so can edit in for loop
  dfname$most_selected_score = 0
  dfname$most_selected_type = 0
  
  #rename columns for nicer legend
  names(dfname)[4:6] = c('missense', 'nonsense', 'splice')
  
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
    
    #adding this label to the row
    dfname$most_selected_type[i] = composite_label
  }
  
  #make sure the label is a factor to allow plotting
  dfname$most_selected_type = as.factor(dfname$most_selected_type)
  
  #create scatter plot, horizontal line indicates the p_value at which we decided the selection was significant
  print(head(dfname))
  
  ggscatter(dfname, x = "most_selected_score", y = "qallsubs_cv", xlab = "Strength of positive selection (for most highly selected mutation type)", ylab = "-log10(P-value)",
            label = "gene_name", repel = TRUE, color = "most_selected_type", size = point_size, font.label = c(font_size, "plain")) + 
    scale_y_continuous(breaks = seq(0, 11, by = 1)) +
    scale_color_discrete(name = 'Most highly selected mutation type(s)') +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red", linewidth = 0.3)
}

#running plots
most_selected_type_plot(slow_genes, 1, 10)
most_selected_type_plot(quick_genes, 0.8, 10)
