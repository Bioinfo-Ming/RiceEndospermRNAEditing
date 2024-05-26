# Author: Ming Chen
# E-mail: chenm@big.ac.cn
# Date: 2023.04.27

# Put libraries here
library(readxl)
library(ggplot2)
library(ggtree)
library(aplot)

# read clustering data
raw_data <- read_excel("E://Rice_Endosperm_RNA_Editing//14.CDS_recoding_synonymous_editing_events_clustering//CDS_recoding_editing_events_clustering.xlsx",
                       sheet = "CDS_recoding_events_clu_all")

head(raw_data)

# give the order of x-axis and y-axis
raw_data$km_model.cluster <- factor(raw_data$km_model.cluster,
                                    levels = c('I', 'II', 'III', 'IV', 'V'),
                                    labels = c('C1', 'C2', 'C3', 'C4', 'C5'))

# plot bubble
cluster_gene_sites <- ggplot(raw_data, aes(x = `km_model.cluster`, y = `Gene`)) +

  geom_point(aes(size = `Standardized_Sites`, color = `Cluster_Number`)) +
  
  scale_y_discrete() +

  theme_bw() +

  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        
        axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14)) +

  # color of bubble
  scale_color_gradient(low = "blue", high = "red") +
  
  # size range of bubble
  scale_size_continuous(range = c(2, 10)) +

  labs(x = "RNA editing sites cluster by K-means",
       y = "") +

  # legend of discrete variable
  guides(size = guide_legend(order = 3))

# data transforming
new_data <- reshape2::dcast(raw_data,
                            Gene~km_model.cluster,
                            value.var = "Standardized_Sites")

rownames(new_data) <- new_data[,1]

new_data <- new_data[,-1]
head(new_data)

row_cluster <- ggtree(hclust(dist(new_data)))

row_cluster + geom_tiplab()



col_data <- reshape2::dcast(raw_data,
                            km_model.cluster~Gene,
                            value.var = "Standardized_Sites")

rownames(col_data) <- col_data[,1]

col_data <- col_data[,-1]

col_cluster <- ggtree(hclust(dist(col_data))) + layout_dendrogram()



cluster_gene_sites%>%insert_left(row_cluster, width = 0.3)
# cluster_gene_sites%>%insert_left(row_cluster, width = 0.3)%>%insert_top(col_cluster, height = 0.15)
# cluster_gene_sites
# cluster_gene_sites%>%insert_left(cluster_gene,width = 0.3)

ggsave(cluster_gene_sites%>%insert_left(row_cluster, width = 0.3),
       filename = 'CDS_recoding_development_clustering_bubble_v4.pdf',
       height = 8,
       width = 6, 
       dpi = 300,
       limitsize = FALSE)

dev.off()