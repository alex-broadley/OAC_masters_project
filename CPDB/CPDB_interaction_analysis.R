## Cell phone DB results analysis

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

#read in pre chemo and post chemo CPDB results
pre_chemo = read.delim("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/CPDB/pre_results/statistical_analysis_significant_means_04_08_2025_102143.txt")
post_chemo = read.delim("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/CPDB/post_results/statistical_analysis_significant_means_04_07_2025_164617.txt")


#function to get dataframe of interactions for a specified cell yupe
extract_interactions = function(df, cellName){
  #extract relevant columns
  general_columns = df[,1:14] 
  cellType_columns = df[, grepl(cellName, names(df)) & !(names(df) %in% names(general_columns))]
  
  #combine relevant columns into dataframe
  combined_df =  cbind(general_columns, cellType_columns)
  
  #identify and remove NA interactions
  interaction_indexes <- apply(combined_df[, 15:ncol(combined_df)], 1, function(row) all(is.na(row)))
  interaction_df = combined_df[interaction_indexes == FALSE, ]
  
  return(interaction_df)
  
}

#get interaction DFs for G0 cells pre and post chemo
G0_pre = extract_interactions(pre_chemo, 'G0')
G0_post = extract_interactions(post_chemo, 'G0')

G0_pre_interactions = unique(G0_pre$id_cp_interaction)
G0_post_interactions = unique(G0_post$id_cp_interaction)

#check overlap between interactions
table(G0_pre_interactions %in% G0_post_interactions)
table(G0_post_interactions %in% G0_pre_interactions)


#get interactions for fast cycling tumo cells
fast_pre = extract_interactions(pre_chemo, 'fast')
fast_post = extract_interactions(post_chemo, 'fast')

fast_pre_interactions = unique(fast_pre$id_cp_interaction)
fast_post_interactions = unique(fast_post$id_cp_interaction)

#check overlap between fast interaction pre/post chemo
table(fast_pre_interactions %in% fast_post_interactions)
table(fast_post_interactions %in% fast_pre_interactions)
