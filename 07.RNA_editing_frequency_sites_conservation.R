# Author: Ming Chen
# E-mail: chenm@big.ac.cn
# Data: 2023.07.17

library(gcookbook)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggrepel)

# load raw data
editing_conservation <- read.table("E:\\Rice_Endosperm_RNA_Editing\\10.RNA_editing_development_evolution_conservation\\11.Rice_Editing_Conservation\\01.species\\RNA_editing_frequency_sites_conservation_species_cyt_c.txt",
                                   header = TRUE)

# plot scatter
scatter_plot <- ggscatter(editing_conservation,
          x = "AverageEditingFrequency",
          y = "EditingAAPercentage") +
  
  # fill color  
  geom_point(aes(color = GeneName), size = 4) +
  
  # fit curve with "glm"
  geom_smooth(method = "glm", color = "black") +
  
  # calculate correlation of "R"
  stat_cor(method = "spearman") +
  
  # add title of x-axis and y-axis
  labs(x = "Editing frequency of CDS-recoding sites",
       y = "Conservation of amino acids",
       size = 20) +
  
  # size of titles
  theme(axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        legend.title = element_blank()) +
  
  theme(panel.grid = element_blank()) +
  
  geom_vline(xintercept = 0.1, lty = "dashed") +
  geom_vline(xintercept = 0.89, lty = "dashed") +
  geom_hline(yintercept = 0.15, lty = "dashed") +
  geom_hline(yintercept = 0.5, lty = "dashed")

# scatter_plot
merged_figure <- scatter_plot + ggrepel::geom_text_repel(aes(label = GeneCDSPosition), editing_conservation)

# merged_figure <- scatter_plot

merged_figure
# ggsave(merged_figure,
#        filename = 'test.pdf',
#        height = 20,
#        width = 40, 
#        dpi = 300,
#        limitsize = FALSE)

# dev.off()