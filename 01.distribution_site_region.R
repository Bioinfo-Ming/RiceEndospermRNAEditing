# chenm@big.ac.cn
# 2023.04.06

library(venn)

venn_data <- read.delim('E:\\Rice_Endosperm_RNA_Editing\\06.distribution_CDS_recoding_C_to_U_editing_sites\\distribution_CDS_recoding_site_region.txt',
                        sep = "\t",
                        header = TRUE,
                        check.names = FALSE)

venn_list <- list(venn_data[,1], venn_data[,2], venn_data[,3], venn_data[,4], venn_data[,5])

names(venn_list) <- colnames(venn_data[1:5])

venn_list <- purrr::map(venn_list, na.omit)

venn(venn_list,
     zcolor = c("#EE3B3B", "#6495ED", "#8B7355", "#EEC900", "#008B8B"),
     opacity = 0.5,
     box = F,
     ilcs = 1,
     sncs = 1.3)