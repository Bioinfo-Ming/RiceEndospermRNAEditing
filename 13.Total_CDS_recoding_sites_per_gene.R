# Author: Ming Chen
# E-mail: chenm@big.ac.cn

# Put libraries here
library(readxl)
library(tidyverse)
library(scales)
library(ggplot2)

# The number of hydrophobic amino acids before and after RNA editing in each gene.
pro_aa_editing <- read_excel("PhoAAEditing_v2.xlsx", sheet = "total")

pro_aa_num_gene <- data.frame(pro_aa_editing$GeneName,
                              pro_aa_editing$StandardizedEditingSites)

show_col(hue_pal()(12))

PhoAABeforeAfterEditing <- ggplot(pro_aa_num_gene, aes(pro_aa_editing$GeneName, pro_aa_editing$StandardizedEditingSites)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.7, fill = "#936eaa") +
  theme_test() +
  
  coord_flip() +
  
  # scale_x_continuous(position = "right") +
  
  theme(axis.text = element_text(colour = 'black', size = 12), 
        legend.title = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank()
        ) +
  
  scale_x_discrete(limits = rev(c('nad1', 'nad2', 'nad4', 'nad4L', 'nad5', 'nad6',
                                'nad7', 'cob', 'cox2', 'atp6', 'atp9', 'ccmB', 
                                'ccmC', 'ccmFc', 'ccmFn', 'rpl2', 'rps1', 'rps13', 
                                'rps19', 'rps4', 'rps7', 'mat-r', 'orf183', 'orf288',
                                'orfX'))) +
  
  labs(x = "Edited Genes", y = "Standardized  CDS-recoding sites per gene")

  
PhoAABeforeAfterEditing
# ggsave('PhoAABeforeAfterEditing.pdf', dpi = 300)