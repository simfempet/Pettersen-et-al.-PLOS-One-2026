library(DESeq2)
library(ggplot2)
library(ggthemes)

Counts <- read.table("C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Submission/Counts.txt", header = T)


rownames(Counts) <- Counts$Geneid

Counts <- Counts[, -c(1:6)]

colnames(Counts)

colnames(Counts) <- c("33C_1", "33C_2", "33C_3",
                      "37C_day1_1", "37C_day1_2", "37C_day1_3",
                      "37C_day7_1", "37C_day7_2", "37C_day7_3",
                      "37C_day14_1", "37C_day14_2", "37C_day14_3",
                      "Ctrl_noFBS_1", "Ctrl_noFBS_2", "Ctrl_noFBS_3",
                      "LPS_1", "LPS_2", "LPS_3",
                      "Ctrl_plusDMSO_1", "Ctrl_plusDMSO_2", "Ctrl_plusDMSO_3",
                      "PP2_1", "PP2_2", "PP2_3")


columns_ordered <- c("33C_1", "33C_2", "33C_3",
                     "37C_day1_1", "37C_day1_2", "37C_day1_3",
                     "37C_day7_1", "37C_day7_2", "37C_day7_3",
                     "37C_day14_1", "37C_day14_2", "37C_day14_3",
                     "Ctrl_noFBS_1", "Ctrl_noFBS_2", "Ctrl_noFBS_3",
                     "Ctrl_plusDMSO_1", "Ctrl_plusDMSO_2", "Ctrl_plusDMSO_3",
                     "LPS_1", "LPS_2", "LPS_3",
                     "PP2_1", "PP2_2", "PP2_3")


Counts <- Counts[, columns_ordered]



#Subset the counts by experiment
Counts_diff <- Counts[, c(1:12)]



Counts_trt <- Counts[, c(13:24)]



#Make coldata for each experiment
colData_diff <- data.frame(label = colnames(Counts_diff),
                      replicate = rep(c("1", "2", "3"), 4),
                      Condition = c(rep("33C", 3),
                                    rep("37C day 1", 3),
                                    rep("37C day 7", 3),
                                    rep("37C day 14", 3)))



colData_trt <- data.frame(label = colnames(Counts_trt),
                      replicate = rep(c("1", "2", "3"), 4),
                      Condition = c(rep("Ctrl - FBS", 3),
                                    rep("Ctrl + DMSO", 3),
                                    rep("LPS trt", 3),
                                    rep("PP2 trt", 3)))



#Make a dds for each experiment
dds_diff <- DESeqDataSetFromMatrix(countData = Counts_diff, colData = colData_diff, design = ~ 1)

keep_diff <- rowSums(counts(dds_diff)) >= 10
dds_diff <- dds_diff[keep_diff, ]
dds_diff <- vst(dds_diff)



dds_trt <- DESeqDataSetFromMatrix(countData = Counts_trt, colData = colData_trt, design = ~ 1)

keep_trt <- rowSums(counts(dds_trt)) >= 10
dds_trt <- dds_trt[keep_trt, ]
dds_trt <- vst(dds_trt)



#Get the counts after filtering and normalization
Counts_diff <- assay(dds_diff)

Counts_trt <- assay(dds_trt)



#PCA
##Diff
PCA_matrix <- Counts_diff



pca <- prcomp(t(PCA_matrix), scale. = T)



pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)



pca.dat <- data.frame(Sample = rownames(pca$x),
                      PC1 = pca$x[,1],
                      PC2 = pca$x[,2],
                      PC3 = pca$x[,3],
                      PC4 = pca$x[,4],
                      Condition = c(rep("33C",3), rep("37C day 1",3), rep("37C day 7", 3), rep("37C day 14", 3)))



pca.dat$Condition <- factor(pca.dat$Condition, levels = c("33C", "37C day 1", "37C day 7", "37C day 14"))



ggplot(pca.dat, aes(x = PC1, y = PC2, fill = Condition))+
  geom_point(pch = 21, size = 3, stroke = 1, color = "black")+
  scale_fill_brewer(palette = "Paired")+
  labs(x = paste("PC1 (", pca.var.per[1], ")%", sep = ""),
       y = paste("PC2 (", pca.var.per[2], ")%", sep = ""))+
  theme_few() + theme(legend.position = c(0.8, 0.3), legend.background = element_rect(colour = "gray"))



#Trt
PCA_matrix <- Counts_trt



pca <- prcomp(t(PCA_matrix), scale. = T)



pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)



pca.dat <- data.frame(Sample = rownames(pca$x),
                      PC1 = pca$x[,1],
                      PC2 = pca$x[,2],
                      PC3 = pca$x[,3],
                      PC4 = pca$x[,4],
                      Condition = c(rep("LPS ctrl",3), rep("PP2 ctrl",3), rep("LPS trt", 3), rep("PP2 trt", 3)))



pca.dat$Condition <- factor(pca.dat$Condition, levels = c("LPS ctrl", "PP2 ctrl", "LPS trt", "PP2 trt"))



ggplot(pca.dat, aes(x = PC1, y = PC2, fill = Condition))+
  geom_point(pch = 21, size = 3, stroke = 1, color = "black")+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[5:8]))+
  labs(x = paste("PC1 (", pca.var.per[1], ")%", sep = ""),
       y = paste("PC2 (", pca.var.per[2], ")%", sep = ""))+
  theme_few() 