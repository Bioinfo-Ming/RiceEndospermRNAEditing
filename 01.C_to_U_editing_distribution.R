# chenm@big.ac.cn
# 2023.03.19

# Average G-to-A editing efficiency during five development stages.
rawdata = read.csv("E:\\Rice_Endosperm_RNA_Editing\\04.analysis_G_to_A_editing\\C_to_U_Editing_Efficience.txt", 
                   header = TRUE, sep = "\t")

# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow = TRUE), height = c(1,9))

# Draw the boxplot and histogram
par(mar = c(0, 4.1, 0.2, 2.1))

boxplot(rawdata$AverageEditing, horizontal = TRUE, ylim = c(0, 1),
        xaxt = "n", frame = F)

par(mar = c(4.5, 4.6, 1.1, 2.1))

hist(rawdata$AverageEditing, breaks = 10, border = F, main = "",
     xlab = "C-to-U average RNA editing efficiency",
     ylab = "Editing sites",
     xlim = c(0, 1),
     ylim = c(0, 70))