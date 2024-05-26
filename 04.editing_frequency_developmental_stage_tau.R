# Author: Ming Chen
# E-mail: chenm@big.ac.cn
# Date: 2023.04.21

library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(rstatix)

# File including RNA editing sites, region and editing frequency at different development stages.
rawdata = read.csv("E:\\Rice_Endosperm_RNA_Editing\\08.distribution_editing_frequency_during_development\\editing_frequency_developmental_stage_tau_v2.txt", 
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

# Order genome regions
rawdata$region <- factor(rawdata$region,
                         levels = c("CDS-recoding", "CDS-synonymous", "Intronic", "Intergenic"))

# Plot basic grouped boxplot.
editing_tau <- ggplot(rawdata, aes(x = region, y = tau)) + 
  
  # add error bar
  stat_boxplot(geom = "errorbar",
               width = 0.2,
               size = 1.3) +
  
  # plot boxplot
  geom_boxplot(aes(fill = region),
               width = 0.4,
               size = 0.8,
               position =  position_dodge(0.75),
               outlier.shape = NA) +
  
  labs(x = "Rice endosperm developmental stages", y = "RNA editing frequency") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  
  # color matching
  ggsci::scale_color_nejm() +

  # filter frey background
  theme_bw()

# add wilcox.test
editing_tau <- editing_tau + geom_signif(comparisons = list(c("CDS-recoding", "CDS-synonymous"),
                                                            c("CDS-recoding", "Intronic"),
                                                            c("CDS-recoding", "Intergenic")),
                                         map_signif_level = T,
                                         test = wilcox.test,
                                         size = 1,
                                         textsize = 6,
                                         y_position = c(1, 1.05, 1.1, 1.15)) +
  
  # graphic beautification
  theme(axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.line = element_line(size = 0.2, colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.03),
        legend.direction = "horizontal",
        legend.text = element_text(size = 12))

editing_tau

# ggsave(editing_site_barplot, filename = 'editing_sites_development_boxplot_v3.pdf', height = 5, width = 8, dpi = 300)