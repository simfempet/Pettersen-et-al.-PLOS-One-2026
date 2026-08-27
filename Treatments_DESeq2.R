library(DESeq2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(ggpubr)
library(ggrepel)
library(ggpmisc)

library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

library(igraph)
library(tidygraph)
library(ggraph)

library(RColorBrewer)

library(ComplexHeatmap)
library(pheatmap)
library(reshape2)

library(GEOquery)
library(GOSemSim)

library(DescTools)


#Organize the count matrix
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

#Subset for treatment samples
Counts_trt <- Counts[, c(13:24)]



#Experimental metadata
colData <- data.frame(label = colnames(Counts_trt),
                      replicate = rep(c("1", "2", "3"), 4),
                      Condition = c(rep("Ctrl - FBS", 3),
                                    rep("Ctrl + DMSO", 3),
                                    rep("LPS trt", 3),
                                    rep("PP2 trt", 3)))

rownames(colData) <- colData$label

all(colnames(Counts_trt) %in% rownames(colData))

all(colnames(Counts_trt) == rownames(colData))



#Construct DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = Counts_trt, colData = colData, design = ~ Condition)



#Pre-filtering: Keep rows that have at least 10 reads in total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]



#set factor level
dds$Condition <- relevel(dds$Condition, ref = "Ctrl - FBS")



#Run DESeq
dds <- DESeq(dds)



#Look at results
#LPS vs ctrl
res_LPS <- results(dds, contrast = c("Condition", "LPS trt", "Ctrl - FBS"), alpha = 0.05)
summary(res_LPS)
res_LPS_df <- data.frame(res_LPS)



#PP2 vs ctrl
res_PP2 <- results(dds, contrast = c("Condition", "PP2 trt", "Ctrl + DMSO"), alpha = 0.05)
summary(res_PP2)
res_PP2_df <- data.frame(res_PP2)



#Mutate the results
res_LPS_df <- res_LPS_df %>% 
  mutate(Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up",
                                log2FoldChange < 0 & padj <= 0.05 ~ "Down",
                                TRUE ~ "Unchanged"),
         GeneID = rownames(res_LPS_df),
         GeneID = sub("\\..*", "", GeneID))



res_PP2_df <- res_PP2_df %>% 
  mutate(Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up",
                                log2FoldChange < 0 & padj <= 0.05 ~ "Down",
                                TRUE ~ "Unchanged"),
         GeneID = rownames(res_PP2_df),
         GeneID = sub("\\..*", "", GeneID))



#Add gene symbols
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res_LPS_df$GeneID,
  mart = ensembl)

res_LPS_df <- merge(res_LPS_df, gene_map, by.x = "GeneID", by.y = "ensembl_gene_id", all.x = TRUE)
res_LPS_df <- res_LPS_df[res_LPS_df$hgnc_symbol != "", ]
res_LPS_df <- res_LPS_df[! is.na(res_LPS_df$Expression), ]



gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res_PP2_df$GeneID,
  mart = ensembl)

res_PP2_df <- merge(res_PP2_df, gene_map, by.x = "GeneID", by.y = "ensembl_gene_id", all.x = TRUE)
res_PP2_df <- res_PP2_df[res_PP2_df$hgnc_symbol != "", ]
res_PP2_df <- res_PP2_df[! is.na(res_PP2_df$Expression), ]



#Get the list of DEGs
LPS__up = (subset(res_LPS_df, res_LPS_df$Expression == "Up"))$hgnc_symbol
LPS__down = (subset(res_LPS_df, res_LPS_df$Expression == "Down"))$hgnc_symbol



PP2__up = (subset(res_PP2_df, res_PP2_df$Expression == "Up" & res_PP2_df$log2FoldChange > 1))$hgnc_symbol
PP2__down = (subset(res_PP2_df, res_PP2_df$Expression == "Down" & res_PP2_df$log2FoldChange < -1))$hgnc_symbol



#GO ORA
##LPS
GO_BP_LPS_up <- enrichGO(gene = LPS__up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_LPS_df$hgnc_symbol)
GO_BP_LPS_up_df <- data.frame(GO_BP_LPS_up)



GO_BP_LPS_down <- enrichGO(gene = LPS__down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_LPS_df$hgnc_symbol)
GO_BP_LPS_down_df <- data.frame(GO_BP_LPS_down)



##PP2
GO_BP_PP2_up <- enrichGO(gene = PP2__up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_PP2_df$hgnc_symbol)
GO_BP_PP2_up_df <- data.frame(GO_BP_PP2_up)



GO_BP_PP2_down <- enrichGO(gene = PP2__down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_PP2_df$hgnc_symbol)
GO_BP_PP2_down_df <- data.frame(GO_BP_PP2_down)



#Volcano plots
label_df_LPS <- res_LPS_df %>%
  filter(Expression != "Unchanged") %>%  
  dplyr::count(Expression) %>%
  mutate(label = paste0("n = ", n), x = ifelse(Expression == "Up", 2, -2), 
         y = max(-log10(ifelse(res_LPS_df$padj == 0, 1e-300, res_LPS_df$padj)), na.rm = TRUE) * 0.7)

ggplot(res_LPS_df, aes(x = log2FoldChange, y = -log(padj, 10)))+
  geom_point(size = 2, alpha = 0.3, pch = 21, aes(fill = Expression))+
  labs(x = "Log2FC ", y = "-Log10FDR", subtitle = "LPS trt vs. Ctrl")+
  scale_fill_manual(values = c("dodgerblue", "gray85", "red3"))+ 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_few()+ ylim(0, 250)+
  geom_label(data = label_df_LPS, aes(x = x, y = y, label = label, fill = Expression),
    color = "black", label.size = 0.5, alpha = 0.5, size = 4, show.legend = F)



label_df_PP2 <- res_PP2_df %>%
  filter(Expression != "Unchanged") %>%  
  dplyr::count(Expression) %>%
  mutate(label = paste0("n = ", n), x = ifelse(Expression == "Up", 2, -2), 
         y = max(-log10(ifelse(res_PP2_df$padj == 0, 1e-300, res_PP2_df$padj)), na.rm = TRUE) * 0.7)

ggplot(res_PP2_df, aes(x = log2FoldChange, y = -log(padj, 10)))+
  geom_point(size = 2, alpha = 0.3, pch = 21, aes(fill = Expression))+
  labs(x = "Log2FC ", y = "-Log10FDR", subtitle = "PP2 trt vs. Ctrl")+
  scale_fill_manual(values = c("dodgerblue", "gray85", "red3"))+ 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_few()+ ylim(0, 250)+
  geom_label(data = label_df_PP2, aes(x = x, y = y, label = label, fill = Expression),
             color = "black", label.size = 0.5, alpha = 0.5, size = 4, show.legend = F)



#Heatmap of selected core genes in PP2 treatment: downregulation of slit diaphragm components and increase in compensation mechanisms
Counts_PP2 <- assay(dds)
Counts_PP2 <- Counts_PP2[,c(4:6, 10:12)]
rownames(Counts_PP2) <- sub("\\..*", "", rownames(Counts_PP2))
Counts_PP2 <- Counts_PP2[res_PP2_df$GeneID, ]



gene_info <- getBM(
  attributes=c("ensembl_gene_id", "transcript_length"),
  filters="ensembl_gene_id",
  values=rownames(Counts_PP2),
  mart=ensembl)



gene_lengths <- aggregate(transcript_length ~ ensembl_gene_id, data=gene_info, median)
rownames(gene_lengths) <- gene_lengths$ensembl_gene_id
common_genes <- intersect(rownames(Counts_PP2), rownames(gene_lengths))
counts <- Counts_PP2[common_genes, ]



lengths <- gene_lengths[common_genes, "transcript_length"]
lengths_kb <- lengths / 1000
rpk <- sweep(counts, 1, lengths_kb, "/")
scaling_factors <- colSums(rpk)
tpm_PP2 <- sweep(rpk, 2, scaling_factors, "/") * 1e6



all(rownames(tpm_PP2) == res_PP2_df$GeneID)
rownames(tpm_PP2) <- res_PP2_df$hgnc_symbol



tpm_PP2_subset <- tpm_PP2[c("MAFB", "TCF21", "FYN", "NCK2",
                            "SYNPO", "CD2AP", "ACTN4", "TJP1", "FAT1", "FAT2"),]



mat <- t(scale(t(tpm_PP2_subset)))



colvector <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "PP2_1", "PP2_2", "PP2_3")



annotation <- data.frame(Condition = c(rep("Ctrl", 3),
                                       rep("PP2 trt", 3)))



colpal <- RColorBrewer::brewer.pal(8, "Paired")[c(6, 8)]



colours <- list('Condition' = c('Ctrl' = colpal[1], 'PP2 trt' = colpal[2]))

colAnn <- HeatmapAnnotation(df = annotation,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),
                            show_legend = T, show_annotation_name = F)



mybreaks <- c(seq(from = min(mat), to = max(mat), by = 0.5))
mycolor <- colorRampPalette(c("dodgerblue", "white", "red3"))(length(mybreaks))



HM_PP2 <- ComplexHeatmap::pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE,
                  color = mycolor, show_rownames = TRUE, 
                  show_colnames = FALSE, labels_col = colvector,
                  name = "TPM Z-score", top_annotation = colAnn)



#Plot LPS results consistent with original study: TRL expression and other cytokines
Counts_LPS <- assay(dds)
Counts_LPS <- Counts_LPS[,c(1:3, 7:9)]
rownames(Counts_LPS) <- sub("\\..*", "", rownames(Counts_LPS))
Counts_LPS <- Counts_LPS[res_LPS_df$GeneID, ]



gene_info <- getBM(
  attributes=c("ensembl_gene_id", "transcript_length"),
  filters="ensembl_gene_id",
  values=rownames(Counts_LPS),
  mart=ensembl)



gene_lengths <- aggregate(transcript_length ~ ensembl_gene_id, data=gene_info, median)
rownames(gene_lengths) <- gene_lengths$ensembl_gene_id
common_genes <- intersect(rownames(Counts_LPS), rownames(gene_lengths))
counts <- Counts_LPS[common_genes, ]



lengths <- gene_lengths[common_genes, "transcript_length"]
lengths_kb <- lengths / 1000
rpk <- sweep(counts, 1, lengths_kb, "/")
scaling_factors <- colSums(rpk)
tpm_LPS <- sweep(rpk, 2, scaling_factors, "/") * 1e6



all(rownames(tpm_LPS) == res_LPS_df$GeneID)
rownames(tpm_LPS) <- res_LPS_df$hgnc_symbol



tpm_LPS_subset1 <- tpm_LPS[c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6"),] #For dotplot

tpm_LPS_subset2 <- tpm_LPS[c("NFKB1", "NFKB2", "RELB", "IL6", "CXCL8", "CCL2", "TNF", 
                             "IL1A", "IL1B", "IL23A", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CCL5", "CSF1", "CSF2", "CSF3"),] #For heatmap



tpm_LPS_subset1 <- melt(tpm_LPS_subset1, varnames = c("Gene", "Sample"), value.name = "TPM")
tpm_LPS_subset1$logTPM <- log2(tpm_LPS_subset1$TPM + 1)
tpm_LPS_subset1$Condition <- tpm_LPS_subset1$Sample
tpm_LPS_subset1$Condition <- sub("_.*", "", tpm_LPS_subset1$Condition)



padj_df <- res_LPS_df[, c("hgnc_symbol", "padj")]
padj_df <- padj_df[padj_df$hgnc_symbol %in% unique(tpm_LPS_subset1$Gene), ]



tpm_LPS_subset1 <- merge(tpm_LPS_subset1, padj_df,
                         by.x = "Gene", by.y = "hgnc_symbol",
                         all.x = TRUE)



padj_positions <- tpm_LPS_subset1 %>%
  group_by(Gene) %>%
  summarise(y_pos = max(logTPM, na.rm = TRUE) + 0.3,
            padj = unique(padj))

padj_positions <- padj_positions %>%
  mutate(stars = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ "ns"))



colpal <- RColorBrewer::brewer.pal(9, "Paired")[5:9]



ggplot(tpm_LPS_subset1, aes(x = Condition, y = logTPM, fill = Condition))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch = 21, alpha = 0.7, dotsize = 1.5)+
  stat_summary(fun = mean, geom = "crossbar", aes(color = Condition))+
  facet_wrap(~ Gene)+
  scale_fill_manual(values = c(colpal[1], colpal[3]))+
  scale_color_manual(values = c(colpal[1], colpal[3]))+
  theme_few() +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none") +
  labs(x = "Condition", y = "log2(TPM + 1)")+
  geom_text(data = padj_positions,
            aes(x = 1.5, y = y_pos, label = stars), 
            inherit.aes = FALSE,
            size = 3.8,
            fontface = "bold")
  


mat2 <- t(scale(t(tpm_LPS_subset2)))



colvector2 <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "LPS_1", "LPS_2", "LPS_3")



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



mybreaks <- c(seq(from = min(mat2), to = max(mat2), by = 0.5))
mycolor <- colorRampPalette(c("dodgerblue", "white", "red3"))(length(mybreaks))



HM_LPS <- ComplexHeatmap::pheatmap(mat2, cluster_cols = TRUE, cluster_rows = TRUE,
                                   color = mycolor, show_rownames = TRUE, 
                                   show_colnames = FALSE, labels_col = colvector2,
                                   name = "TPM Z-score", top_annotation = colAnn)



#Podocyte barrier score
Counts_VST <- assay(dds)
Counts_VST <- varianceStabilizingTransformation(Counts_VST)
rownames(Counts_VST) <- sub("\\..*", "", rownames(Counts_VST))
Counts_VST <- Counts_VST[res_PP2_df$GeneID, ]



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



VST_sub <- Counts_VST[c("KIRREL1", "SYNPO", "CD2AP", "FAT1", "FAT2", "TJP1", "TJP2",
                        "TRPC6", "CDH3"),]



mat <- t(scale(t(VST_sub)))



colnames(mat)[1:6] <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4", "Ctrl_5", "Ctrl_6")



mybreaks <- c(seq(from = min(mat), to = max(mat), by = 0.5))
mycolor <- colorRampPalette(c("dodgerblue", "white", "red3"))(length(mybreaks))



HM <- ComplexHeatmap::pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE,
                               color = mycolor, show_rownames = TRUE, 
                               show_colnames = TRUE, name = "Z-score")



score_df <- data.frame(Sample = colnames(mat),
                       Score = colMeans(mat),
                       Condition = c(rep("Ctrl", 6),
                                     rep("LPS trt", 3),
                                     rep("PP2 trt", 3)))



score_df$Sample <- factor(score_df$Sample,
                          levels = c("PP2_1", "PP2_2", "PP2_3",
                                     "Ctrl_3", "Ctrl_1", "Ctrl_4", "Ctrl_6", "Ctrl_2", "Ctrl_5",
                                     "LPS_1", "LPS_2", "LPS_3"))



score_lm <- lm(Score ~ Condition, score_df)
summary(score_lm)
DunnettTest(Score ~ Condition, score_df, control = "Ctrl")



ggplot(score_df, aes(x = Condition, y = Score, fill = Condition))+
  geom_dotplot(binaxis='y', stackdir='center', pch=21, color="black", dotsize=1.5)+
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(5,7,8)])+
  scale_color_manual(values=brewer.pal(12, "Paired")[c(5,7,8)])+
  stat_summary(fun=mean, geom="crossbar", aes(color = Condition))+
  stat_compare_means(comparisons=list(c("Ctrl","LPS trt"), c("Ctrl", "PP2 trt")), method="t.test", label="p.signif")+
  theme_few()+ theme(legend.position="none")+
  labs(x = NULL, y = "Podocyte barrier score")+ ylim(-1.1, 1.1)



#Analysis of hallmark gene sets related to unspecific toxicity response
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



ranks_PP2 <- res_PP2_df$stat
names(ranks_PP2) <- res_PP2_df$hgnc_symbol
ranks_PP2 <- ranks_PP2[!is.na(names(ranks_PP2))]
ranks_PP2 <- sort(ranks_PP2, decreasing = TRUE)



fgsea_PP2 <- fgsea(pathways = pathway_list, stats = ranks_PP2, minSize = 15, maxSize = 500)
fgsea_PP2 <- fgsea_PP2 %>%
  mutate(Significance = case_when(NES > 0 & padj <= 0.05 ~ "Up",
                                  NES < 0 & padj <= 0.05 ~ "Down",
                                  TRUE ~ "Unchanged"))



src_pathway_only <- pathway_list["GAUTSCHI_SRC_SIGNALING"]
fgsea_src_direct <- fgsea(pathways = src_pathway_only, stats = ranks_PP2, minSize = 8, maxSize = 500)



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



ranks_LPS <- res_LPS_df$stat
names(ranks_LPS) <- res_LPS_df$hgnc_symbol
ranks_LPS <- ranks_LPS[!is.na(names(ranks_LPS))]
ranks_LPS <- sort(ranks_LPS, decreasing = TRUE)



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



heat_df <- bind_rows(fgsea_PP2 %>%
    filter(pathway %in% stress_terms) %>%
    transmute(pathway, NES, padj, Treatment = "PP2 (75 µM)"),
    fgsea_LPS %>%
    filter(pathway %in% stress_terms) %>%
    transmute(pathway, NES, padj, Treatment = "LPS (40 µg/mL)")) %>%
  mutate(sig = case_when(padj < 0.001 ~ "***", 
                         padj < 0.01  ~ "**",
                         padj < 0.05  ~ "*",
                         TRUE ~ ""))



heat_df$pathway <- gsub("HALLMARK_", "", heat_df$pathway)
heat_df$pathway <- gsub("GOBP_", "", heat_df$pathway)
heat_df$pathway <- gsub("_", " ", heat_df$pathway)



pathway_order <- heat_df %>%
  filter(Treatment == "PP2 (75 µM)") %>%
  arrange(desc(NES)) %>%
  pull(pathway)

heat_df$pathway <- factor(heat_df$pathway, levels = pathway_order)



ggplot(heat_df, aes(x = pathway, y = Treatment, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 5, fontface = "bold") +
  scale_fill_distiller(palette = "RdBu", name = "NES\n(Trt vs. Ctrl)", limits = c(-3, 3.5))+
  labs(x = NULL, y = NULL) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        legend.position = "top")



#Add on: Western blot quantification heatmap
library(readxl)
WB_quant <- read_excel("Skole/UiB/MCD-model/RNA-seq/WB_quant.xlsx", 
                       sheet = "Sheet4")

WB_quant$Condition <- factor(WB_quant$Condition, levels = c("Ctrl", "PP2 trt"))
WB_quant$Lane <- factor(WB_quant$Lane, levels = c("NPHS1", "pNPHS1"))

ggplot(WB_quant, aes(x = Condition, y = Lane, fill = Normalized_Int))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "#AA3377", midpoint = 0.5, name = "Normalized band intensity")+
  theme_few()+ theme(text = element_blank())



#Save outputs (for semantic similarity analysis later)
saveRDS(res_LPS_df, file = "RNA-seq/DE_LPS.RDS")
saveRDS(res_PP2_df, file = "RNA-seq/DE_PP2.RDS")
saveRDS(GO_BP_LPS_up, file = "RNA-seq/GO_LPS_up.RDS")
saveRDS(GO_BP_LPS_down, file = "RNA-seq/GO_LPS_down.RDS")
saveRDS(GO_BP_PP2_up, file = "RNA-seq/GO_PP2_up.RDS")
saveRDS(GO_BP_PP2_down, file = "RNA-seq/GO_PP2_down.RDS")