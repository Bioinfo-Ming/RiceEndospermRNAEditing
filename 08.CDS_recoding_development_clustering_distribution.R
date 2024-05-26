# Author: Ming Chen
# E-mail: chenm@big.ac.cn

library(ggplot2)
library(ggpubr)

# input data
editing_frequency_distribution <- read.table("E:\\Rice_Endosperm_RNA_Editing\\14.CDS_recoding_editing_events_clustering\\CDS_recoding_development_clustering_distribution.txt",
                                             header = TRUE, sep = "\t")

# sort info
editing_frequency_distribution$Cluster_N <- factor(editing_frequency_distribution$Cluster_N,
                                                   levels = c('4', '1', '2', '3', '5'))

editing_frequency_distribution$Developmental_Stage <- factor(editing_frequency_distribution$Developmental_Stage,
                                                             levels = c('3_DAF', '6_DAF', '9_DAF',
                                                                        '12_DAF', '15_DAF'))
head(editing_frequency_distribution)

# plot boxplot and line
editing_frequency_distribution <- ggplot(editing_frequency_distribution,
                                         aes(x = Developmental_Stage,
                                             y = Editing_Frequency)) +
   
  geom_line(aes(group = Genome_Position), color = "grey", size = 0.5) +
  
  # add error bar
  stat_boxplot(geom = "errorbar",
               width = 0.5,
               size = 0.5,
               color = "#4393C3") +
  
  # add x labels
  scale_x_discrete(labels = c("3", "6", "9", "12", "15")) +
  
  # add x + y titles
  labs(x = "Developmental stages", y = "RNA editing frequency", size = 16) +
  
  # plot boxplot
  geom_boxplot(aes(fill = Cluster_N), fill = "white", color = "#4393C3", size = 0.5) +
  
  # add lines of mean values
  stat_summary(fun = mean, geom = "line", aes(group = Cluster_N), color = "red",
               size = 0.8) +
  
  # add points of mean values
  stat_summary(fun = mean, geom = "point", color = "red", shape = 17, fill = "red", 
               size = 3) +
  
  # facet
  facet_wrap(vars(Cluster_N), nrow = 1) +
  
  # delete grey background
  theme_bw() +
  
  theme(
        panel.grid = element_blank(),
        
        # delete background of title
        strip.background = element_blank(),
        
        # delete title of each facet
        strip.text.x = element_blank()) +
  
  # change scale of y-axis
  scale_y_continuous(name = "RNA editing frequency",
                     breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))

editing_frequency_distribution
# ggsave("CDS_recoding_development_clustering_distribution.pdf", dpi = 300)
# dev.off()