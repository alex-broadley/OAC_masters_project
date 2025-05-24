setwd('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data')

library(ggplot2)
library(pBrackets)
library(stringr)

### Stacked bar plot for pooled pre- and post-chemo proportion of tumour cell types ------------------------

#read counts data in from scanpy
t_plot_data = read.csv('G0_plot_agreement_quartiles.csv')
length(t_plot_data)
t_plot_data$donor_chemo = paste(t_plot_data$donor, t_plot_data$chemo)
t_plot_data$chemo = factor(t_plot_data$chemo, levels = c('pre', 'post'))

#create stacked bar plot - label by count and proportion on each bar
pre_post_plot = ggplot(t_plot_data, aes(x = chemo, y = X._by_treatment, fill = G0_class))+
  geom_col(color = "black")+
  scale_color_manual(values=c("#c2c2d6", "#e60039", "#3366ff"))+
  scale_fill_manual(values=c("#c2c2d6", "#e60039", "#3366ff")) +
  geom_text(size = 5, fontface = 'bold', 
            aes(label = paste('   ', round(X._by_treatment, 2),"%", '\n n =', nb_by_treatment)), 
            position = position_stack(vjust = 0.5))+
  labs(
    x = "Chemotherapy stage",
    y = "% of total tumour cells by chemo stage"
  ) +
  theme_classic() +
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position = 'none') + #remove legened as will be shared with next plot
  ylim(0, 110) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))

### fishers exact test ------ G0 pre vs post--------------------------------------------------------------

#counts of G0 and non-G0 cells before and after chemo
dat <- data.frame(
  "pre_chemo" = c(212, 626),
  "post_chemo" = c(16, 58),
  row.names = c("G0", "Not-G0"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("pre_chemo", "post_chemo")

#this is >5, which tutorial said means there is sufficient statistical power to use fishers
chisq.test(dat)$expected

#fishers exact test comparing the two, p = 0.5758 (< 0.05) so is not significant
fishersExact <- fisher.test(dat)

### fishers exact test ------ fast-cycling pre vs post----------------------------------------------------

#counts of fast and non-fast cells before and after chemo
dat <- data.frame(
  "pre_chemo" = c(202, 636),
  "post_chemo" = c(26, 48),
  row.names = c("G0", "Not-G0"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("pre_chemo", "post_chemo")

#this is >5, which tutorial said means there is sufficient statistical power to use fishers
chisq.test(dat)$expected

#fishers exact test comparing the two, p = 0.02773 (< 0.05) so is significant
fishersExact <- fisher.test(dat, alternative = 'less')

### Stacked bar plot for pooled pre- and post-chemo proportion of tumour cell types ------------------------

#read in data from scanpy
patient_plot_data = read.csv('G0_patient_plotData.csv', row.names = 'X')

#change prolif labels to capitals for better legend
patient_plot_data$G0_class = str_to_sentence(patient_plot_data$G0_class)

library(dplyr)

#label sample by chemo stage
patient_plot_data <- patient_plot_data %>%
  mutate(chemo = case_when(
    donor %in% c("Pt7", "Pt5", "Pt4", "Pt2") ~ "Pre-chemo",
    donor %in% c("Pt8", "Pt1", "Pt6", "Pt3") ~ "Post-chemo"
  ))

#order samples by their chemo stage so they are plotted in order: pre-chemo then post-chemo
patient_plot_data <- patient_plot_data %>%
  mutate(donor = factor(donor, levels = unique(donor[order(chemo, decreasing = TRUE)])))

#create stacked bar plot - same labelling as above
patient_plot = ggplot(patient_plot_data, aes(x = donor, y = X._by_patient, fill = G0_class))+
  geom_col(color = "black")+
  scale_color_manual(values=c("#c2c2d6", "#e60039", "#3366ff"))+
  scale_fill_manual(values=c("#c2c2d6", "#e60039", "#3366ff"), name = 'Proliferation type') +
  geom_text(size = 5, fontface = 'bold', 
            aes(label = ifelse(nb_by_patient == 0, NA, 
                               paste('   ', round(X._by_patient, 2), "%", '\n n =', nb_by_patient))),
            position = position_stack(vjust = 0.5)
  ) +
  labs(
    x = "Sample",
    y = "% of total tumour cells by patient"
  ) +
  theme_classic() +
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size = 15)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))
  

library(gridExtra)

# put the two figures next to each other, make plot by patient wider to accomodate extra bars
figure = gridExtrafigure = grid.arrange(pre_post_plot, patient_plot, ncol = 2, widths = c(0.6, 2))



