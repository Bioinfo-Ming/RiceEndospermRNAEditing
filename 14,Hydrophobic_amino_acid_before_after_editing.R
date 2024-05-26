# Author: Ming Chen
# E-mail: chenm@big.ac.cn

# Put libraries here
library(readxl)
library(tidyverse)
library(scales)
library(ggplot2)

# The number of hydrophobic amino acids before and after RNA editing in each gene.
pro_aa_editing <- read_excel("PhoAAEditing_v2.xlsx", sheet = "pho")

pro_aa_num_gene <- data.frame(pro_aa_editing$GeneName, pro_aa_editing$EditingState, 
                              pro_aa_editing$StandardizedEditingSites)

show_col(hue_pal()(12))

ggplot(pro_aa_num_gene, aes(pro_aa_editing$GeneName, pro_aa_editing$StandardizedEditingSites, 
                            fill = pro_aa_editing$EditingState)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +
  theme_test() +
  
  theme(axis.text = element_text(colour = 'black', size = 12), 
        legend.title = element_blank(), 
        panel.grid = element_blank(), 
        ) +
  
  scale_x_discrete(limits = rev(c('nad1', 'nad2', 'nad4', 'nad4L', 'nad5', 'nad6',
                                'nad7', 'cob', 'cox2', 'atp6', 'atp9', 'ccmB', 
                                'ccmC', 'ccmFc', 'ccmFn', 'rpl2', 'rps1', 'rps13', 
                                'rps19', 'rps4', 'rps7', 'mat-r', 'orf183', 'orf288',
                                'orfX'))) +
  
  labs(x = "Edited Genes", y = "Hydrophobic amino acids before and after RNA editing") +

  coord_flip()

ggsave('PhoAABeforeAfterEditing.pdf', dpi = 300)