# Author: Ming Chen
# E-mail: chenm@big.ac.cn
# Date: 2023.06.27

library(ggplot2)

# file including RNA editing types, developmental stages and RNA editing sites
rawdata <- read.table("E:\\Rice_Endosperm_RNA_Editing\\18.Editing_sites_distribution_at_5_developmental_stages\\all_editing_sites_at_5_developmental_stages.txt", 
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

# order RNA editing types
rawdata$EditingType <- factor(rawdata$EditingType,
                              levels = c('C-to-U (MT)',
                                         'G-to-A (MT)',
                                         'C-to-U (PT)',
                                         'U-to-C (MT)'))

rawdata$DevelopmentalStage <- factor(rawdata$DevelopmentalStage,
                                     levels = c('3_DAF',
                                                '6_DAF',
                                                '9_DAF',
                                                '12_DAF',
                                                '15_DAF'))

clustered_bars <- ggplot(data = rawdata,
                         aes(x = DevelopmentalStage,
                             y = EditingSites,
                             fill = EditingType)) +
  
  geom_bar(stat = "identity",
           position = "dodge",
           width = 0.9) +
  
  labs(x = "Rice endosperm development stages",
       y = "RNA editing sites") +
  
  # scale_fill_manual(values = c("#BCE6D8", "#F8F8F8", "#DFE7EA", "#4D5174")) +
  
  # delete color of background
  theme_classic() +
  
  theme(axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = c(0.93, 0.88),
        legend.text = element_text(size = 10))

clustered_bars