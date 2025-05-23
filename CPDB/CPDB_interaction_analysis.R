## Cell phone DB results analysis

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")


pre_chemo1 = read.delim("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/CPDB/pre_results/statistical_analysis_significant_means_04_08_2025_102143.txt")
pre_chemo2 = read.delim("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/CPDB/pre_results/statistical_analysis_significant_means_04_08_2025_105627.txt")

pre_chemo = rbind(pre_chemo1, pre_chemo2)


post_chemo1 = read.delim("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/CPDB/post_results/statistical_analysis_significant_means_04_07_2025_164617.txt")
post_chemo2 = read.delim("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/CPDB/post_results/statistical_analysis_significant_means_04_08_2025_094429.txt")

post_chemo = rbind(post_chemo1, post_chemo2)


extract_interactions = function(df, cellName){
  general_columns = df[,1:14] 
  cellType_columns = df[, grepl(cellName, names(df)) & !(names(df) %in% names(general_columns))]
  combined_df =  cbind(general_columns, cellType_columns)
  
  interaction_indexes <- apply(combined_df[, 15:ncol(combined_df)], 1, function(row) all(is.na(row)))
  
  interaction_df = combined_df[interaction_indexes == FALSE, ]
  
  return(interaction_df)
  
}

G0_pre = extract_interactions(pre_chemo, 'G0')
G0_post = extract_interactions(post_chemo, 'G0')

G0_pre_interactions = unique(G0_pre$id_cp_interaction)
G0_post_interactions = unique(G0_post$id_cp_interaction)

table(G0_pre_interactions %in% G0_post_interactions)
table(G0_post_interactions %in% G0_pre_interactions)


fast_pre = extract_interactions(pre_chemo, 'fast')
fast_post = extract_interactions(post_chemo, 'fast')

fast_pre_interactions = unique(fast_pre$id_cp_interaction)
fast_post_interactions = unique(fast_post$id_cp_interaction)

table(fast_pre_interactions %in% fast_post_interactions)
table(fast_post_interactions %in% fast_pre_interactions)
