# Survival analysis based on bulk RNA classification of slow or aggressive tumors
library(ggpubr)
library(survminer)
library(survival)

setwd("/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data")

#adapted from Maria's code 


# load bulk file w/ attached G0 classifications so can group patients by proliferation score
bulk_RNA = read.csv('bulk_RNA_classified.csv', row.names = 'X')

#load in survival data
survival_data = read.csv('fullset_occams_SA_SAZ_28052020.csv', row.names = 'X')

mapping = bulk_RNA[, c("DNAid", "proliferation")]

survival_classified = merge(survival_data, mapping, by.x = "Illumina.ID", by.y = 'DNAid')
#fit with standard function and cox proportional hazard model
fit <- survfit( Surv(Weeks.Survival.c) ~ proliferation, data = survival_classified )
fit2 <- coxph( Surv(Weeks.Survival.c) ~ proliferation, data = survival_classified )

# Survival plot -> directly from Maria code but changed color
ggsurvplot(
  fit,
  data = survival_classified,
  size = 1,                 # change line size
  palette =
    c("#e60039","#3366ff"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when
  ggtheme = theme_bw()      # Change ggplot2 theme
)

#alternative survival curve plot
ggadjustedcurves(fit2, data = survival_classified, 
                 variable = "proliferation",
                 individual.curves=TRUE)

