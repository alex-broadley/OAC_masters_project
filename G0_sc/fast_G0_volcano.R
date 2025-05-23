library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)


fast_G0_DE = read.csv('/Users/alex/Documents/BIOL0041-Project/OAC_masters_project/data/fast_G0_DE.csv')

fast_G0_DE$sig_direction = "Not signifcant"

fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges > 0.5 & fast_G0_DE$pvals_adj < 0.05] <- "Upregulated"
fast_G0_DE$sig_direction[fast_G0_DE$logfoldchanges < -0.5 & fast_G0_DE$pvals_adj < 0.05] <- "Downregulated"

fast_G0_DE$delabel = NA

for (row in 1:nrow(fast_G0_DE)) {
  if (fast_G0_DE$names[row] %in% head(fast_G0_DE[order(fast_G0_DE$pvals_adj), "names"], 50)){
    print(fast_G0_DE$names[row])
    fast_G0_DE$delabel[row] = fast_G0_DE$names[row]
  }
  else {
    fast_G0_DE$delabel[row] = NA
  }
}

#### create volcano plot 
ggplot(data = fast_G0_DE, aes(x = logfoldchanges, y = -log10(pvals_adj), col = sig_direction, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +
  geom_point() +
  scale_color_manual(values = c("darkblue", "grey", "red"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_classic() +
  geom_text_repel(max_overlaps = "50") + 
  scale_fill_discrete(name = "")
