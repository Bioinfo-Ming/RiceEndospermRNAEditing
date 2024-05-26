# Author: Ming Chen
# E-mail: chenm@big.ac.cn
# Date: 2023.04.21

library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(rstatix)

# file including RNA editing sites, region and editing frequency at different development stages.
rawdata <- read.table("E:\\Rice_Endosperm_RNA_Editing\\08.distribution_editing_frequency_during_development\\editing_frequency_developmental_stage_v2.txt", 
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

# convert wide data into long data.
data_long <- pivot_longer(rawdata,
                         cols =! region,
                         names_to = "developmental_stage",
                         values_to = "editing_frequency")

# sort information
data_long$developmental_stage <- factor(data_long$developmental_stage,
                                        levels = c('3_DAF', '6_DAF', '9_DAF', '12_DAF', '15_DAF'),
                                        labels = c('3', '6', '9', '12', '15'))

data_long$region <- factor(data_long$region,
                           levels = c('CDS-recoding', 'CDS-synonymous', 'Intronic', 'Intergenic'))

# Plot basic grouped boxplot.
editing_site_barplot <- ggplot(data_long,
                               aes(x = developmental_stage,
                                   y = editing_frequency)) +
  # add error bar
  stat_boxplot(geom = "errorbar",
               width = 0.25,
               size = 0.8,
               aes(fill = region),
               position = position_dodge(0.75)) +
  
  # divide into groups
  geom_boxplot(aes(fill = region),
               width = 0.6,
               position =  position_dodge(0.75),
               
               # delete outlier
               outlier.shape = NA) +

  # # add points of mean values
  # stat_summary(fun = mean,
  #              geom = "point",
  #              color = "red",
  #              shape = 17,
  #              fill = "red",
  #              size = 2,
  #              aes(group = region),
  #              position = position_dodge(0.75)) +

  # add titles
  labs(x = "Rice endosperm developmental stages",
       y = "RNA editing frequency") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +

  # color matching
  ggsci::scale_color_nejm() +

  # filter frey background
  theme_bw()

# wilcox_test
stat.test <- data_long %>%
  group_by(developmental_stage) %>%
  wilcox_test(
    editing_frequency ~ region,
    comparisons = list(c("CDS-recoding", "CDS-synonymous"),
                       c("CDS-recoding", "Intronic"),
                       c("CDS-recoding", "Intergenic")),
    paired = FALSE,
    alternative = "two.sided",
    p.adjust.method = "bonferroni")

# Determine the position (level) of p-value
stat.test <- stat.test %>% add_xy_position(x = "developmental_stage")

editing_site_barplot <- editing_site_barplot +
  
  stat_pvalue_manual(stat.test,
                     y.position = c(1.03, 1.07, 1.11),
                     label = "p.adj.signif",
                     size = 5,
                     
                     # width of line
                     bracket.size = 1,
                     
                     # length of line
                     tip.length = 0.01) +

  # graphic beautification
  theme(axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.line = element_line(size = 0.2, colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.5, 1),
        legend.direction = "horizontal",
        legend.text = element_text(size = 12))

editing_site_barplot
# ggsave(editing_site_barplot, filename = 'editing_frequency_developmental_stage.pdf', height = 5, width = 8, dpi = 300)