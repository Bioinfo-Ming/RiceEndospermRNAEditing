# Author: Ming Chen
# E-mail: chenm@big.ac.cn
# Date: 2023.08.29

library(ggplot2)
library(viridis)

# Input raw data.
editing_events <- read.table("E:\\Rice_Endosperm_RNA_Editing\\19.CDS_recoding_synonymous_editing_events_during_development\\CDS_recoding_synonymous_intronic_overlappping_editing_events_development.txt",
                             header = TRUE)

# Plot line chart.
ggplot(editing_events, 
       mapping = aes(x = `DevelopmentStage`,
                     y = `EditingFrequency`,
                     color = `Region`,
                     group = `Position`,
                     shape = `Gene`)) +
  
  geom_line() +
  
  facet_wrap(vars(`Gene`), nrow = 5, strip.position = "top") +
  
  # Delete the grey background.
  theme_bw() +
  
  scale_x_continuous(breaks = seq(0, 15, 3)) +
  
  theme(panel.grid = element_blank(),
        legend.position="right") +
  
  xlab("Rice endosperm developmental stages") +
  ylab("RNA editing frequency")