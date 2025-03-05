setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

library(ggplot2)

### Stacked bar plot for pooled pre- and post-chemo proportion of tumour cell types

t_plot_data = read.csv('G0_treatment_plotData.csv')
length(t_plot_data)
t_plot_data$donor_chemo = paste(t_plot_data$donor, t_plot_data$chemo)

colnames(t_plot_data)[3] = 'Cycle_type'

t_plot_data$chemo = factor(t_plot_data$chemo, levels = c('pre', 'post'))



ggplot(t_plot_data, aes(x = chemo, y = X._by_treatment, fill = Cycle_type))+
  geom_col()+
  scale_color_manual(values=c("#286837", "#772989", "#bbbbbb"))+
  scale_fill_manual(values=c("#286837", "#772989", "#bbbbbb")) +
  geom_text(aes(label = paste('   ', round(X._by_treatment, 2),"%", '\n n =', nb_by_treatment)), position = position_stack(vjust = 0.5))+
  labs(
    x = "Chemotherapy stage",
    y = "% of total tumour cells"
  ) +
  theme_classic()

### fishers exact test

#counts of G0 and non-G0 cells before and after chemo
dat <- data.frame(
  "pre_chemo" = c(59, 197),
  "post_chemo" = c(31, 52),
  row.names = c("G0", "Not-G0"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("pre_chemo", "post_chemo")

#this is >5, which tutorial said means there is sufficient statistical power to use fishers
chisq.test(dat)$expected

#fishers exact test comparing the two, p = 0.01463 (< 0.05) so is significant
fishersExact <- fisher.test(dat)

#Chi-squared test 
X_sq = chisq.test(dat)

