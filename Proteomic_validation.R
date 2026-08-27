library(DEP)
library(dplyr)
library(tidyverse)
library(reshape2)
library(stringr)
library(DESeq2)
library(DescTools)

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggExtra)

library(ComplexHeatmap)

library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)

library(RColorBrewer)

library(fgsea)
library(decoupleR)
library(biomaRt)



colpal <- c("#FB9A99", "#FDBF6F", "#FF7F00")



#Load data
proteinGroups <- read.table("C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Submission/proteinGroups.txt", header = T)
data <- proteinGroups



LFQ_columns <- grep("LFQ.", colnames(data))
colnames(data)[LFQ_columns]



#Filter out contaminants
data <- subset(data, data$Potential.contaminant != "+") 



#Remove duplicated gene names
data$Gene.names %>% duplicated() %>% any()



duplicates <- data %>% 
  group_by(Gene.names) %>%
  summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% 
  filter(frequency > 1)



data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")



data_unique$name %>% duplicated() %>% any()


#Experimental metadata
Experimental_design <- data.frame(condition = c(rep("Ctrl", 3),
                                                rep("LPS_trt", 3),
                                                rep("PP2_trt", 3)),
                                  replicate = c(1:3, 1:3, 1:3),
                                  label = gsub("^LFQ\\.intensity\\.", "", colnames(data)[LFQ_columns]))



rownames(Experimental_design) <- Experimental_design$label



LFQ_columns <- grep("LFQ.", colnames(data_unique)) 



LFQ_matrix <- data_unique[, LFQ_columns] 



colnames(LFQ_matrix) = str_remove(colnames(LFQ_matrix), "LFQ.intensity.") 
all(colnames(LFQ_matrix) %in% rownames(Experimental_design))
all(colnames(LFQ_matrix) == rownames(Experimental_design))



data_se <- make_se(data_unique, LFQ_columns, Experimental_design)



#Filter proteins: thresh = 3/3 samples within at least one condition (Force NPHS1 to be included)
exprs_mat <- assay(data_se)



conds <- colData(data_se)$condition



keep_protein <- apply(exprs_mat, 1, function(x) {
  tapply(!is.na(x), conds, function(v) mean(v) == (1))
})

keep_rows <- apply(keep_protein, 2, any)
keep_rows["NPHS1"] <- TRUE



data_filt <- data_se[keep_rows, ]



protein_counts <- colSums(!is.na(assay(data_filt)))



df_counts <- data.frame(Sample = names(protein_counts),
                        Proteins = protein_counts,
                        Condition = colData(data_filt)$condition)
df_counts$Condition <- factor(df_counts$Condition)



ggplot(df_counts, aes(x = Sample, y = Proteins, fill = Condition))+
  geom_bar(stat = "identity", color = "black")+theme_few()+ 
  scale_fill_manual(values = colpal, name = "Condition")+
  labs(x = NULL, y = "Number of proteins")+
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")



#Normalize the data
data_norm <- normalize_vsn(data_filt)



meanSdPlot(data_norm)$gg + theme_few() + scale_fill_viridis_c(option = "rocket")+ labs(x = "Rank", y = "SD")



#Impute missing values
plot_detect(data_filt)



set.seed(2156)
data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)



norm_df <- as.data.frame(assay(data_norm))
norm_df_long <- norm_df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
norm_df_long$assay <- "Normalized"



imp_df <- as.data.frame(assay(data_imp))
imp_df_long <- imp_df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
imp_df_long$assay <- "Imputed"



long_df <- rbind(norm_df_long, imp_df_long)



ggplot(long_df, aes(x = value, group = assay, fill = assay))+
  geom_density(alpha = 0.5, color = "black")+
  scale_fill_manual(values = c("gray35", "gray95"), name = "Layer")+
  labs(x = "Log2 LFQ intensity", y = "Density")+
  theme_few()



#RLI plot
LFQ_matrix_imp <- as.matrix(assay(data_imp))



RLI_matrix <- LFQ_matrix_imp



RLI_matrix <- sweep(RLI_matrix, MARGIN = 1, STATS = rowMedians(RLI_matrix))



RLI_df <- melt(RLI_matrix, varnames = c("Protein", "Sample"), value.name = "RLI")



RLI_df$Condition <- RLI_df$Sample
RLI_df$Condition <- ifelse(grepl("^Ctrl_[^_]+$", RLI_df$Condition), 
                           "Ctrl", 
                           sub("^([^_]+_[^_]+)_.*$", "\\1", RLI_df$Condition))
RLI_df$Condition <- factor(RLI_df$Condition)



ggplot(data = RLI_df, aes(x = Sample, y = RLI, fill = Condition))+
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2)+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = colpal)+
  labs(x = NULL, y = "Relative log intensity")+
  theme_few() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")



#PCA
gene_var <- apply(LFQ_matrix_imp, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:1000]



PCA_matrix <- LFQ_matrix_imp[top_genes,]



pca <- prcomp(t(PCA_matrix), scale. = T)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)



pca.dat <- data.frame(Sample = rownames(pca$x),
                      X = pca$x[,1],
                      Y = pca$x[,2])
pca.dat$condition <- Experimental_design$condition
pca.dat$condition <- factor(pca.dat$condition)



ggplot(pca.dat, aes(x = X, y = Y, fill = condition))+
  geom_point(pch = 21, size = 3, stroke = 1, color = "black")+
  scale_fill_manual(values = colpal, name = "Condition", labels = c("Ctrl", "LPS trt", "PP2 trt"))+
  labs(x = paste("PC1 (", pca.var.per[1], ")%", sep = ""),
       y = paste("PC2 (", pca.var.per[2], ")%", sep = ""))+
  theme_few()



#DEPs
data_diff_PP2 <- test_diff(data_imp, type = "manual", test = "PP2_trt_vs_Ctrl")
dep_PP2 <- add_rejections(data_diff_PP2, alpha = 0.05, lfc = 0)
PP2_results <- get_results(dep_PP2)

PP2_results <- PP2_results %>%
  mutate(Expression = case_when(significant == "TRUE" & PP2_trt_vs_Ctrl_ratio > 0 ~ "Up",
                                significant == "TRUE" & PP2_trt_vs_Ctrl_ratio < 0 ~ "Down",
                                TRUE ~ "Unchanged"))



data_diff_LPS <- test_diff(data_imp, type = "manual", test = "LPS_trt_vs_Ctrl")
dep_LPS <- add_rejections(data_diff_LPS, alpha = 0.05, lfc = 0)
LPS_results <- get_results(dep_LPS)

LPS_results <- LPS_results %>%
  mutate(Expression = case_when(significant == "TRUE" & LPS_trt_vs_Ctrl_ratio > 0 ~ "Up",
                                significant == "TRUE" & LPS_trt_vs_Ctrl_ratio < 0 ~ "Down",
                                TRUE ~ "Unchanged"))



#Volcano plots
label_df <- PP2_results %>%
  filter(Expression != "Unchanged") %>%  
  dplyr::count(Expression) %>%
  mutate(label = paste0("n = ", n), x = ifelse(Expression == "Up", 2, -2), 
         y = max(-log10(ifelse(PP2_results$PP2_trt_vs_Ctrl_p.adj == 0, 1e-300, PP2_results$PP2_trt_vs_Ctrl_p.adj)), na.rm = TRUE) * 0.7)



p1 <- ggplot(PP2_results, aes(x = PP2_trt_vs_Ctrl_ratio, y = -log(PP2_trt_vs_Ctrl_p.adj, 10)))+
  geom_point(size = 2, alpha = 0.3, pch = 21, aes(fill = Expression))+
  labs(x = "Log2FC ", y = "-Log10FDR", subtitle = "PP2 trt vs. Ctrl")+
  scale_fill_manual(values = c("dodgerblue", "gray85", "red3"))+ 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_few()+ theme(legend.position = "none")+
  geom_label(data = label_df, aes(x = x, y = y, label = label, fill = Expression),
             color = "black", label.size = 0.5, alpha = 0.5, size = 4, show.legend = F, fontface = "bold")



label_df <- LPS_results %>%
  filter(Expression != "Unchanged") %>%  
  dplyr::count(Expression) %>%
  mutate(label = paste0("n = ", n), x = ifelse(Expression == "Up", 2, -2), 
         y = max(-log10(ifelse(LPS_results$LPS_trt_vs_Ctrl_p.adj == 0, 1e-300, LPS_results$LPS_trt_vs_Ctrl_p.adj)), na.rm = TRUE) * 0.7)



p2 <- ggplot(LPS_results, aes(x = LPS_trt_vs_Ctrl_ratio, y = -log(LPS_trt_vs_Ctrl_p.adj, 10)))+
  geom_point(size = 2, alpha = 0.3, pch = 21, aes(fill = Expression))+
  labs(x = "Log2FC ", y = "-Log10FDR", subtitle = "LPS trt vs. Ctrl")+
  scale_fill_manual(values = c("dodgerblue", "gray85", "red3"))+ 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_few()+ theme(legend.position = "none")+
  geom_label(data = label_df, aes(x = x, y = y, label = label, fill = Expression),
             color = "black", label.size = 0.5, alpha = 0.5, size = 4, show.legend = F, fontface = "bold")



ggarrange(p1, p2)



#Get the list of DEGs
LPS__up = (subset(LPS_results, LPS_results$Expression == "Up"))$name
LPS__down = (subset(LPS_results, LPS_results$Expression == "Down"))$name



PP2__up = (subset(PP2_results, PP2_results$Expression == "Up"))$name
PP2__down = (subset(PP2_results, PP2_results$Expression == "Down"))$name



#GO ORA
##LPS
GO_BP_LPS_up <- enrichGO(gene = LPS__up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = LPS_results$name)
GO_BP_LPS_up_df <- data.frame(GO_BP_LPS_up)



GO_BP_LPS_down <- enrichGO(gene = LPS__down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = LPS_results$name)
GO_BP_LPS_down_df <- data.frame(GO_BP_LPS_down)



##PP2
GO_BP_PP2_up <- enrichGO(gene = PP2__up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = PP2_results$name)
GO_BP_PP2_up_df <- data.frame(GO_BP_PP2_up)



GO_BP_PP2_down <- enrichGO(gene = PP2__down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = PP2_results$name)
GO_BP_PP2_down_df <- data.frame(GO_BP_PP2_down)



#Look at NPHS1 and NPHS2 (PP2 vs ctrl)
GOIs <- c("NPHS1", "NPHS2")



LFQ_matrix_sub <- LFQ_matrix_imp[GOIs, c(1:3, 7:9)]



LFQ_df <- as.data.frame(LFQ_matrix_sub)



long_df <- LFQ_df %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "label",        
    values_to = "Expression")
long_df$Condition <- rep(c(rep("Ctrl", 3), rep("PP2 trt", 3)), 2)
long_df$Condition <- factor(long_df$Condition)



NPHS1_lm <- lm(Expression ~ Condition,
               subset(long_df, long_df$Gene == "NPHS1"))
summary(NPHS1_lm)



NPHS2_lm <- lm(Expression ~ Condition,
               subset(long_df, long_df$Gene == "NPHS2"))
summary(NPHS2_lm)



colpal2 <- c("#E31A1C", "#FF7F00")
ggplot(long_df, aes(x = Gene, y = Expression, fill = Condition))+
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(0.5),
               pch = 21, dotsize = 1.5, stackgroups = TRUE, alpha = 0.5)+
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.5),
               aes(group = Condition, color = Condition), width = 0.4)+
  stat_compare_means(aes(group = Condition), method = "t.test", label = "p.signif")+
  scale_fill_manual(values = colpal2)+
  scale_color_manual(values = colpal2)+
  theme_few()+ ylim(20, 27)+
  labs(x = NULL, y = "Log2 intensity")



#Inflammation markers
GOIs <- c("TNFAIP3", "SAA2", "PYCARD", "ALPK1", "RGS10", "NRROS", "CRADD")



mat <- t(scale(t(LFQ_matrix_imp[GOIs, c(1:6)])))



colvector <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "LPS_1", "LPS_2", "LPS_3")



annotation <- data.frame(Condition = c(rep("Ctrl", 3),
                                       rep("LPS trt", 3)))



colpal <- RColorBrewer::brewer.pal(8, "Paired")[c(5, 7)]

colours <- list('Condition' = c('Ctrl' = colpal[1], 'LPS trt' = colpal[2]))



colAnn <- HeatmapAnnotation(df = annotation,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),
                            show_legend = T, show_annotation_name = F)



mybreaks <- c(seq(from = min(mat), to = max(mat), by = 0.5))
mycolor <- colorRampPalette(c("dodgerblue", "white", "red3"))(length(mybreaks))



HM_LPS <- ComplexHeatmap::pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE,
                                   color = mycolor, show_rownames = TRUE, 
                                   show_colnames = FALSE, labels_col = colvector,
                                   name = "Z-score", top_annotation = colAnn)



#PP2 GOIs
GOIs <- c("FYN", "NCK2", "SYNPO", "CD2AP", "FAT1", "TJP1")



LFQ_matrix_sub <- LFQ_matrix_imp[GOIs, c(1:3, 7:9)]



mat <- t(scale(t(LFQ_matrix_sub)))



LFQ_df <- as.data.frame(mat)



long_df <- LFQ_df %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "label",        
    values_to = "Expression")
long_df$Condition <- rep(c(rep("Ctrl", 3), rep("PP2 trt", 3)), 6)
long_df$Condition <- factor(long_df$Condition)
long_df$Gene <- factor(long_df$Gene, levels = rev(c("FYN", "NCK2", "SYNPO", "CD2AP", "FAT1", "TJP1")))



p1 <- ggplot(long_df, aes(x = Gene, y = Expression, fill = Condition))+
  stat_summary(fun.data = mean_sd, geom = "errorbar", position = position_dodge(0.5),
               aes(group = Condition, color = Condition), width = 0.4)+
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(0.5),
               aes(group = Condition, color = Condition), width = 0.4)+
  scale_fill_manual(values = colpal2)+
  scale_color_manual(values = colpal2)+
  theme_few()+ theme(legend.position = "top")+
  labs(x = NULL, y = "Scaled intensity")+ coord_flip()



diff_df <- long_df %>%
  group_by(Gene, Condition) %>%
  summarize(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = Mean_Expression) %>%
  mutate(Difference = `PP2 trt` - Ctrl)
diff_df$Gene <- factor(diff_df$Gene, levels = rev(c("FYN", "NCK2", "SYNPO", "CD2AP", "FAT1", "TJP1")))



p2 <- ggplot(diff_df, aes(y= Gene, x = " ", fill = Difference))+
  geom_tile(color = "gray35")+
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red3",
                      midpoint = 0, name = "Mean difference\n(PP2 - Ctrl)")+
  theme_void()+ theme(legend.position = "top")



ggarrange(p1, p2, ncol = 2, align = "hv", widths = c(1, 0.7))



#GSEA for unspecific toxcicity response
library(fgsea)
library(msigdbr)



hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")



celldeath_df <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_id %in% c("M46836", "M47138", "M24370", "M46862"))



msigdb_extra_df <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_id %in% c("M46836", "M47138", "M24370", "M46862") | 
           gs_name == "GAUTSCHI_SRC_SIGNALING")



combined_df <- bind_rows(hallmark_df, msigdb_extra_df)



pathway_list <- combined_df %>%
  split(x = .$gene_symbol, f = .$gs_name)



ranks_PP2 <- PP2_results$PP2_trt_vs_Ctrl_ratio
names(ranks_PP2) <- PP2_results$name
ranks_PP2 <- ranks_PP2[!is.na(names(ranks_PP2))]
ranks_PP2 <- sort(ranks_PP2, decreasing = T)



fgsea_PP2 <- fgsea(pathways = pathway_list, stats = ranks_PP2, minSize = 15, maxSize = 500)
fgsea_PP2 <- fgsea_PP2 %>%
  mutate(Significance = case_when(NES > 0 & padj <= 0.05 ~ "Up",
                                  NES < 0 & padj <= 0.05 ~ "Down",
                                  TRUE ~ "Unchanged"))



src_pathway_only <- pathway_list["GAUTSCHI_SRC_SIGNALING"]
fgsea_src_direct <- fgsea(pathways = src_pathway_only, stats = ranks_PP2, minSize = 1, maxSize = 500)



stress_terms <- c(
  "HALLMARK_APOPTOSIS",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_DNA_REPAIR",
  "GOBP_FERROPTOSIS",
  "GOBP_NECROPTOTIC_PROCESS",
  "GOBP_PYROPTOTIC_INFLAMMATORY_RESPONSE")



fgsea_PP2_sub <- fgsea_PP2 %>%
  filter(pathway %in% stress_terms) %>%
  arrange(padj)



fgsea_PP2_sub <- fgsea_PP2_sub[order(fgsea_PP2_sub$padj), ]
fgsea_PP2_sub$pathway <- factor(fgsea_PP2_sub$pathway, levels = fgsea_PP2_sub$pathway)



ggplot(fgsea_PP2_sub, aes(x = -log10(padj), y = pathway, color = NES))+
  geom_point(size = 4)+
  geom_segment(aes(x = 0, xend = -log10(padj), y = pathway, yend = pathway))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray55")+
  scale_color_distiller(palette = "RdBu")+
  labs(x = "-log10(FDR)", y = NULL, color = "NES\n(PP2 (75 µM) vs Ctrl)")+
  theme_few()+ theme(legend.position = "top")


ranks_LPS <- LPS_results$LPS_trt_vs_Ctrl_ratio
names(ranks_LPS) <- LPS_results$name
ranks_LPS <- ranks_LPS[!is.na(names(ranks_LPS))]
ranks_LPS <- sort(ranks_LPS, decreasing = T)



fgsea_LPS <- fgsea(pathways = pathway_list, stats = ranks_LPS, minSize = 15, maxSize = 500)
fgsea_LPS <- fgsea_LPS %>%
  mutate(Significance = case_when(NES > 0 & padj <= 0.05 ~ "Up",
                                  NES < 0 & padj <= 0.05 ~ "Down",
                                  TRUE ~ "Unchanged"))



fgsea_LPS_sub <- fgsea_LPS %>%
  filter(pathway %in% stress_terms) %>%
  arrange(padj)



fgsea_LPS_sub <- fgsea_LPS_sub[order(fgsea_LPS_sub$padj), ]
fgsea_LPS_sub$pathway <- factor(fgsea_LPS_sub$pathway, levels = fgsea_LPS_sub$pathway)



ggplot(fgsea_LPS_sub, aes(x = -log10(padj), y = pathway, color = NES))+
  geom_point(size = 4)+
  geom_segment(aes(x = 0, xend = -log10(padj), y = pathway, yend = pathway))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray55")+
  scale_color_distiller(palette = "RdBu")+
  labs(x = "-log10(FDR)", y = NULL, color = "NES\n(LPS (40 µg/mL) vs Ctrl)")+
  theme_few() + theme(legend.position = "top")



#Barrier score
POIs <- c("KIRREL1", "SYNPO", "CD2AP", "FAT1", "FAT2", "TJP1", "TJP2","TRPC6", "CDH3")
which(POIs %in% rownames(LFQ_matrix_imp))
POIs <- POIs[c(2, 3, 4, 5, 6, 7)]



LFQ_sub <- LFQ_matrix_imp[POIs,]



mat <- t(scale(t(LFQ_sub)))



score_df <- data.frame(Sample = colnames(mat),
                       Score = colMeans(mat),
                       Condition = Experimental_design$condition)
score_df$Condition <- factor(score_df$Condition)



summary(lm(Score ~ Condition, score_df))



ggplot(score_df, aes(x = Condition, y = Score, fill = Condition))+
  geom_dotplot(binaxis='y', stackdir='center', pch=21, color="black", dotsize=1.5)+
  stat_summary(fun=mean, geom="crossbar", aes(color = Condition))+
  theme_few()+ theme(legend.position="none")+
  labs(x = NULL, y = "Podocyte barrier score")+
  stat_compare_means(comparisons = list(c("Ctrl", "LPS_trt"), c("Ctrl", "PP2_trt")), method = "t.test", label = "p.signif")+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+ ylim(-1, 1)+
  scale_x_discrete(labels = c("Ctrl", "LPS trt", "PP2 trt"))



#TF activity inference
library(effsize)
library(effectsize)



net <- read.csv("C:/Users/SimenFemangerPetters/Downloads/CollecTRI.csv")



TF_activity <- run_ulm(mat = LFQ_matrix_imp[, c(1:3, 7:9)], network = net, .source = "source", .target = "target", .mor = "weight")



TF_activity <- subset(TF_activity, TF_activity$source == "TCF21")



TF_mat <- TF_activity %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()



TF_mat <- scale(TF_mat)



plot_df <- data.frame(Sample = rownames(TF_mat),
                      Score = TF_mat[, 1],
                      Condition = c(rep("Ctrl", 3), rep("PP2 trt", 3)))



ggplot(plot_df, aes(x = Condition, y = Score, fill = Condition))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch = 21, dotsize = 1.5)+
  stat_summary(fun = mean, geom = "crossbar", aes(color = Condition))+
  scale_fill_manual(values = colpal2)+
  scale_color_manual(values = colpal2)+
  theme_few()+ theme(legend.position = "none")+
  labs(x = NULL, y = "TCF21 activity")+
  annotate("text", label = "Hedges' g: 1.40 [-0.23, 2.93]", x = 1.5, y = 2)
  


summary(lm(Score ~ Condition, plot_df))
cohen.d(Score ~ Condition, plot_df)
hedges_g(plot_df$Score[4:6], plot_df$Score[1:3])



#Overall alignment between transcriptome and proteome
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
Counts <- Counts[, c(13:24)]
colData <- data.frame(label = colnames(Counts),
                      replicate = rep(c("1", "2", "3"), 4),
                      Condition = c(rep("Ctrl - FBS", 3),
                                    rep("Ctrl + DMSO", 3),
                                    rep("LPS trt", 3),
                                    rep("PP2 trt", 3)))
rownames(colData) <- colData$label
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = colData, design = ~ Condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
Counts_VST <- varianceStabilizingTransformation(assay(dds))
rownames(Counts_VST) <- sub("\\..*", "", rownames(Counts_VST))



ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(Counts_VST),
  mart = ensembl)
gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
Counts_VST <- Counts_VST[rownames(Counts_VST) %in% gene_map$ensembl_gene_id, ]
gene_map <- gene_map[match(rownames(Counts_VST), gene_map$ensembl_gene_id), ]
rownames(Counts_VST) <- gene_map$hgnc_symbol

common_genes <- intersect(rownames(LFQ_matrix_imp), rownames(Counts_VST))

LFQ_matrix_imp <- LFQ_matrix_imp[common_genes, ]
Counts_VST <- Counts_VST[common_genes, ]

all(rownames(LFQ_matrix_imp) == rownames(Counts_VST))



#Plot
library(ggpmisc)
plot_df <- bind_rows(
  data.frame(Gene = names(rowMeans(Counts_VST[, 1:6])), 
             Transcriptome = rowMeans(Counts_VST[, 1:6]), 
             Proteome = rowMeans(LFQ_matrix_imp[, 1:3]), 
             Condition = "Ctrl"),
  data.frame(Gene = names(rowMeans(Counts_VST[, 7:9])), 
             Transcriptome = rowMeans(Counts_VST[, 7:9]), 
             Proteome = rowMeans(LFQ_matrix_imp[, 4:6]), 
             Condition = "LPS"),
  data.frame(Gene = names(rowMeans(Counts_VST[, 10:12])), 
             Transcriptome = rowMeans(Counts_VST[, 10:12]), 
             Proteome = rowMeans(LFQ_matrix_imp[, 7:9]), 
             Condition = "PP2")
)



ggplot(plot_df, aes(x = Transcriptome, y = Proteome, color = Condition))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, color = "gray30")+
  stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE, color = "gray30")+
  scale_color_manual(values = colpal)+
  facet_wrap(~ Condition)+
  labs(x = "Transcriptome\n(normalized counts)", y = "Proteome\n(Log2 intensity)")+
  theme_few()+ theme(legend.position = "none")