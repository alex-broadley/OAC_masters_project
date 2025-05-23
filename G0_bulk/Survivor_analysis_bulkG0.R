# Survival analysis based on bulk RNA classification of slow or aggressive tumors
library(ggpubr)
library(survminer)
library(survival)

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

# load bulk file w/ attached G0 classifications -> take both the raw counts file and the normalised as classifications in normalised
bulk_RNA = read.csv('bulk_RNA_classified.csv', row.names = 'X')

survival_data = read.csv('fullset_occams_SA_SAZ_28052020.csv', row.names = 'X')

mapping = bulk_RNA[, c("DNAid", "proliferation")]

survival_classified = merge(survival_data, mapping, by.x = "Illumina.ID", by.y = 'DNAid')

fit <- survfit( Surv(Weeks.Survival.c) ~ proliferation, data = survival_classified )

fit2 <- coxph( Surv(Weeks.Survival.c) ~ proliferation, data = survival_classified )

ggsurvplot(
  fit,
  data = survival_classified,
  size = 1,                 # change line size
  palette =
    c("#BE92A2","#575D90"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)


ggadjustedcurves(fit2, data = survival_classified, 
                 variable = "proliferation",
                 individual.curves=TRUE)

