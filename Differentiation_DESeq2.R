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
Counts <- Counts_genocode_v26_UiB.S03.2025

rownames(Counts) <- Counts$Geneid

Counts <- Counts[, -c(1:6)]

colnames(Counts)

colnames(Counts) <- c("33C_1", "33C_2", "33C_3",
                      "37C_day1_1", "37C_day1_2", "37C_day1_3",
                      "37C_day14_1", "37C_day14_2", "37C_day14_3",
                      "37C_day7_1", "37C_day7_2", "37C_day7_3",
                      "Chip_1", "Chip_2", "Chip_3", "Chip_4", "Chip_5", "Chip_6",
                      "Ctrl_noFBS_1", "Ctrl_noFBS_2", "Ctrl_noFBS_3",
                      "Ctrl_plusDMSO_1", "Ctrl_plusDMSO_2", "Ctrl_plusDMSO_3",
                      "LPS_1", "LPS_2", "LPS_3",
                      "PP2_1", "PP2_2", "PP2_3")


columns_ordered <- c("33C_1", "33C_2", "33C_3",
                    "37C_day1_1", "37C_day1_2", "37C_day1_3",
                    "37C_day7_1", "37C_day7_2", "37C_day7_3",
                    "37C_day14_1", "37C_day14_2", "37C_day14_3",
                    "Ctrl_noFBS_1", "Ctrl_noFBS_2", "Ctrl_noFBS_3",
                    "Ctrl_plusDMSO_1", "Ctrl_plusDMSO_2", "Ctrl_plusDMSO_3",
                    "LPS_1", "LPS_2", "LPS_3",
                    "PP2_1", "PP2_2", "PP2_3",
                    "Chip_1", "Chip_2", "Chip_3", "Chip_4", "Chip_5", "Chip_6")

Counts <- Counts[, columns_ordered]


#Subset for the differentiation samples
Counts_diff <- Counts[, c(1:12)]



#Experimental metadata
colData <- data.frame(label = colnames(Counts_diff),
                      replicate = rep(c("1", "2", "3"), 4),
                      Condition = c(rep("33C", 3),
                                    rep("37C day 1", 3),
                                    rep("37C day 7", 3),
                                    rep("37C day 14", 3)))
rownames(colData) <- colData$label

all(colnames(Counts_diff) %in% rownames(colData))

all(colnames(Counts_diff) == rownames(colData))



#Construct DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = Counts_diff, colData = colData, design = ~ Condition)



#Pre-filtering: Keep rows that have at least 10 reads in total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]



#set factor level
dds$Condition <- relevel(dds$Condition, ref = "37C day 14")



#Run DESeq
dds <- DESeq(dds)


#Fetch results
res_UndiffvsDiff <- results(dds, contrast = c("Condition", "37C day 14", "33C"), alpha = 0.05)

summary(res_UndiffvsDiff)



#Plot volcano
res_UndiffvsDiff_df <- as.data.frame(res_UndiffvsDiff)
res_UndiffvsDiff_df$GeneID <- rownames(res_UndiffvsDiff_df)
res_UndiffvsDiff_df$GeneID <- sub("\\..*", "", res_UndiffvsDiff_df$GeneID)



ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res_UndiffvsDiff_df$GeneID,
  mart = ensembl)

res_UndiffvsDiff_df <- merge(res_UndiffvsDiff_df, gene_map, by.x = "GeneID", by.y = "ensembl_gene_id", all.x = TRUE)
res_UndiffvsDiff_df <- res_UndiffvsDiff_df[res_UndiffvsDiff_df$hgnc_symbol != "", ]



res_UndiffvsDiff_df <- res_UndiffvsDiff_df %>% 
  mutate(Expression = case_when(log2FoldChange > 0.5 & padj <= 0.05 ~ "Up",
                                log2FoldChange < -0.5 & padj <= 0.05 ~ "Down",
                                TRUE ~ "Unchanged"))



label_df <- res_UndiffvsDiff_df %>%
  filter(Expression != "Unchanged") %>%  
  dplyr::count(Expression) %>%
  mutate(label = paste0("n = ", n), x = ifelse(Expression == "Up", 2, -2), 
         y = max(-log10(ifelse(res_UndiffvsDiff_df$padj == 0, 1e-300, res_UndiffvsDiff_df$padj)), na.rm = TRUE) * 0.7)



ggplot(res_UndiffvsDiff_df, aes(x = log2FoldChange, y = -log(padj, 10)))+
  geom_point(size = 2, alpha = 0.3, pch = 21, aes(fill = Expression))+
  labs(x = "Log2FC ", y = "-Log10FDR", subtitle = "37C day 14 vs. 33C")+
  scale_fill_manual(values = c("dodgerblue", "gray85", "red3"))+ 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(0.5, -0.5), linetype = "dashed")+
  theme_few()+ theme(legend.position = c(0.88, 0.75), legend.background = element_rect(colour = "gray"))+
  geom_label(data = label_df, aes(x = x, y = y, label = label, fill = Expression),
             color = "black", label.size = 0.5, alpha = 0.5, size = 4, show.legend = F, fontface = "bold")



nrow(subset(res_UndiffvsDiff_df, res_UndiffvsDiff_df$Expression != "Unchanged")) /nrow(res_UndiffvsDiff_df) #Fraction of DEGs



#Get gene sets of interest
#Podocyte markers
cell_markers <- read.delim(file = "C:/Users/SimenFemangerPetters/Downloads/PanglaoDB_markers_27_Mar_2020.tsv.gz")
podocyte_markers <- subset(cell_markers, cell_markers$cell.type == "Podocytes" & cell_markers$canonical.marker == "1")

#DEGs
Sigs_up <- (subset(res_UndiffvsDiff_df, res_UndiffvsDiff_df$Expression == "Up"))$hgnc_symbol

Sigs_down <- (subset(res_UndiffvsDiff_df, res_UndiffvsDiff_df$Expression == "Down"))$hgnc_symbol



#Gene Ontology ORA
GO_BP_Up <- enrichGO(gene = Sigs_up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_UndiffvsDiff_df$hgnc_symbol)
GO_BP_Up_df <- as.data.frame(GO_BP_Up)

GO_BP_Down <- enrichGO(gene = Sigs_down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_UndiffvsDiff_df$hgnc_symbol)
GO_BP_Down_df <- as.data.frame(GO_BP_Down)



#Plotting GO results
#Step 1: subset out key terms
keywords <- c("glomerulus", "glomerular", "kidney", "podocyte") #Upregulated terms
pattern <- paste(keywords, collapse = "|")



GO_BP_Up_df <- GO_BP_Up_df[grepl(pattern, GO_BP_Up_df$Description, ignore.case = TRUE), ]
GO_BP_Up_df <- GO_BP_Up_df[order(GO_BP_Up_df$FoldEnrichment, decreasing = T),]
GO_BP_Up_df$Description <- factor(GO_BP_Up_df$Description, levels = rev(GO_BP_Up_df$Description))



ggplot(GO_BP_Up_df, aes(x = FoldEnrichment, y = Description, fill = FoldEnrichment))+
  geom_bar(stat = "identity", color = "black")+  
  scale_fill_distiller(palette = "Reds", direction = 1)+
  geom_label(aes(label = Description),
            hjust = ifelse(GO_BP_Up_df$FoldEnrichment > 0, 1.02, -0.02),
            fill = "white", color = "black", size = 3, fontface = "bold", label.size = 0.3)+
  labs(x = "Fold enrichment" ,y = NULL)+
  theme_few()+ theme(legend.position = "none", axis.text.y = element_blank())



GO_BP_Up@result <- GO_BP_Up@result[grepl(pattern, GO_BP_Up@result$Description, ignore.case = TRUE), ]
GO_BP_Up <- pairwise_termsim(GO_BP_Up)



keywords <- c("cell cycle", "cell division") #Downregulated terms
pattern <- paste(keywords, collapse = "|")



GO_BP_Down_df <- GO_BP_Down_df[grepl(pattern, GO_BP_Down_df$Description, ignore.case = TRUE), ]
GO_BP_Down_df <- GO_BP_Down_df[order(GO_BP_Down_df$FoldEnrichment, decreasing = T),]
GO_BP_Down_df$Description <- factor(GO_BP_Down_df$Description, levels = rev(GO_BP_Down_df$Description))



ggplot(GO_BP_Down_df[c(1, 6, 7, 11, 13), ], aes(x = FoldEnrichment, y = Description, fill = FoldEnrichment))+
  geom_bar(stat = "identity", color = "black")+  
  scale_fill_distiller(palette = "Blues", direction = 1)+
  geom_label(aes(label = Description,
                 hjust = ifelse(FoldEnrichment > 0, 1.02, -0.02)),  
             fill = "white", color = "black", size = 3, fontface = "bold", label.size = 0.3) +
  labs(x = "Fold enrichment" ,y = NULL)+
  theme_few()+ theme(legend.position = "none", axis.text.y = element_blank())



GO_BP_Down@result <- GO_BP_Down@result[grepl(pattern, GO_BP_Down@result$Description, ignore.case = TRUE), ]
GO_BP_Down <- pairwise_termsim(GO_BP_Down)
GO_BP_Down@result$FoldEnrichment <- GO_BP_Down@result$FoldEnrichment * -1



#Step 2: Create enrichment map
GO_combined <- GO_BP_Up
GO_combined@result <- bind_rows(GO_BP_Up@result, GO_BP_Down@result)



keywords <- c("negative")
pattern <- paste(keywords, collapse = "|")



GO_combined@result <- GO_combined@result[!grepl(pattern, GO_combined@result$Description, ignore.case = TRUE), ]



GO_combined <- pairwise_termsim(GO_combined)



nodes <- as.data.frame(GO_combined)

edges <- as.data.frame(as.table(GO_combined@termsim)) %>%
  filter(!is.na(Freq), Freq > 0.2, Var1 != Var2) %>%
  rename(from = Var1, to = Var2, weight = Freq)

nodes <- nodes %>%
  rename(name = Description)

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)



ggraph(graph, layout = "fr") +
  geom_edge_link(aes(alpha = weight), show.legend = FALSE, color = "gray45") +
  geom_node_point(aes(size = Count, color = FoldEnrichment)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_gradient2(name = "Relative fold enrichment\n(37C day 14 vs 33C)",
                        low = "dodgerblue", mid = "white", high = "red3")+
  theme_void()



#Add a clustering algorithm 
graph <- graph %>% 
  mutate(module = as.factor(group_louvain()))

module_labels <- graph %>% 
  as_tibble() %>% 
  group_by(module) %>% 
  summarise(label = name[which.max(Count)])

module_labels[3,2] <- "kidney development\nglomerulus development\npodocyte differentiation"

graph <- graph %>% 
  left_join(module_labels, by = "module")



ggraph(graph, layout = "kk") +
  geom_edge_link(aes(alpha = weight), show.legend = FALSE, color = "gray45") +
  geom_node_point(aes(size = Count, color = FoldEnrichment)) +
  geom_node_text(data = . %>% as_tibble() %>% group_by(module) %>% slice_max(Count, n = 1),
    aes(label = label), size = 4, fontface = "bold", color = "black", repel = TRUE, bg.color = "white")+
  scale_color_gradient2(name = "Relative fold enrichment\n(37C day 14 vs 33C)",
                        low = "dodgerblue", mid = "white", high = "red3")+
  theme_void()



#Heatmap of the podocyte markers of interest
markers <- intersect(podocyte_markers$official.gene.symbol, Sigs_up)



Counts <- assay(dds)
rownames(Counts) <- sub("\\..*", "", rownames(Counts))



res_UndiffvsDiff_df <- res_UndiffvsDiff_df[! is.na(res_UndiffvsDiff_df$GeneID), ]
Counts <- Counts[res_UndiffvsDiff_df$GeneID, ]



ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


gene_info <- getBM(
  attributes=c("ensembl_gene_id", "transcript_length"),
  filters="ensembl_gene_id",
  values=rownames(Counts),
  mart=ensembl)



gene_lengths <- aggregate(transcript_length ~ ensembl_gene_id, data=gene_info, median)
rownames(gene_lengths) <- gene_lengths$ensembl_gene_id

common_genes <- intersect(rownames(Counts), rownames(gene_lengths))
counts <- Counts[common_genes, ]



lengths <- gene_lengths[common_genes, "transcript_length"]
lengths_kb <- lengths / 1000
rpk <- sweep(counts, 1, lengths_kb, "/")
scaling_factors <- colSums(rpk)

tpm <- sweep(rpk, 2, scaling_factors, "/") * 1e6



gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(tpm),
  mart = ensembl)

gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]

tpm <- tpm[rownames(tpm) %in% gene_map$ensembl_gene_id, ]

gene_map <- gene_map[match(rownames(tpm), gene_map$ensembl_gene_id), ]

rownames(tpm) <- gene_map$hgnc_symbol



tpm_subset <- tpm[markers, ]



mat <- t(scale(t(tpm_subset)))

mybreaks <- c(seq(from = min(mat), to = max(mat), by = 0.1))
mycolor <- colorRampPalette(c("dodgerblue", "white", "red3"))(length(mybreaks))



annotation <- data.frame(Condition = c(rep("33C", 3),
                                       rep("37C day 1", 3),
                                       rep("37C day 7", 3),
                                       rep("37C day 14", 3)))

colpal <- RColorBrewer::brewer.pal(4, "Paired")

colours <- list('Condition' = c('33C' = colpal[1], '37C day 1' = colpal[2], 
                                '37C day 7' = colpal[3], '37C day 14' = colpal[4]))

colAnn <- HeatmapAnnotation(df = annotation,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),
                            show_legend = T, show_annotation_name = F)



annotation2 <- data.frame(Protein = rownames(tpm_subset))
annotation2$Podocyte_marker <- 1

Master_TF <- c("MAFB", "FOXC1")
Slit_diaphragm <- c("SYNPO", "MAGI2", "RAB3B", "PLEKHH2", "EPB41L5")
GBM_component <- c("LAMB2", "LAMA3", "COL4A4", "COL4A5")
Cell_adhesion <-  c("ITGA3", "ITGB5", "CDH13")
ECM_modulation <- c("PLOD2", "LOXL2")
Actin_related <- c("MYO1E", "KANK1", "ARHGAP24", "RHPN1", "AIF1L", "CLIC5")
Signaling <- c("VEGFA", "ANGPTL2", "ADM", "PROS1", "MERTK", "NPR3", "PTH1R", "ADCY1")
Stress_protection <- c("OPTN", "LCN2", "HSPB8")

annotation2 <- annotation2 %>%
  mutate(Slit_diaphragm = case_when(Protein %in% Slit_diaphragm ~ 1, TRUE ~ 0),
         Master_TF = case_when(Protein %in% Master_TF ~ 1, TRUE ~ 0),
         GBM_component = case_when(Protein %in% GBM_component ~ 1, TRUE ~ 0),
         Cell_adhesion = case_when(Protein %in% Cell_adhesion ~ 1, TRUE ~ 0),
         ECM_modulation = case_when(Protein %in% ECM_modulation ~ 1, TRUE ~ 0),
         Actin_related = case_when(Protein %in% Actin_related ~ 1, TRUE ~ 0),
         Signaling = case_when(Protein %in% Signaling ~ 1, TRUE ~ 0),
         Stress_protection = case_when(Protein %in% Stress_protection ~ 1, TRUE ~ 0))

annotation2_mat <- annotation2[, -1]
rownames(annotation2_mat) <- annotation2$Protein
annotation2_mat <- annotation2_mat[rownames(mat), ]

rowAnn <- rowAnnotation(
  df = annotation2_mat,
  col = list(
    Podocyte_marker = c("0" = "white", "1" = "black"),
    Master_TF = c("0" = "white", "1" = "black"),
    Slit_diaphragm = c("0" = "white", "1" = "black"),
    GBM_component = c("0" = "white", "1" = "black"),
    Cell_adhesion = c("0" = "white", "1" = "black"),
    ECM_modulation = c("0" = "white", "1" = "black"),
    Actin_related = c("0" = "white", "1" = "black"),
    Signaling = c("0" = "white", "1" = "black"),
    Stress_protection = c("0" = "white", "1" = "black")
  ),
  show_legend = FALSE,
  show_annotation_name = TRUE,      
  annotation_name_rot = 75,         
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 7.5, fontface = "bold"))



HM <- ComplexHeatmap::pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE,
                              show_rownames = TRUE, show_colnames = FALSE,
                              name = "TPM Z-score", color = mycolor,
                              top_annotation = colAnn, left_annotation = rowAnn,
                              fontsize_row = 7.5)



#TF activity inference
library(decoupleR)
net <- read.csv("C:/Users/SimenFemangerPetters/Downloads/CollecTRI.csv")



Counts_VST <- assay(dds)
Counts_VST <- varianceStabilizingTransformation(Counts_VST)
rownames(Counts_VST) <- sub("\\..*", "", rownames(Counts_VST))
Counts_VST <- Counts_VST[res_UndiffvsDiff_df$GeneID, ]



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



TF_activity <- run_ulm(mat = Counts_VST, network = net, .source = "source", .target = "target", .mor = "weight")

TF_activity_sub <- TF_activity %>%
  filter(source %in% c("WT1", "MAFB", "TCF21", "FOXC1", "FOXC2")) %>% #Podocyte master TFs
  mutate(Time = case_when(condition %in% c("37C_day14_1", "37C_day14_2", "37C_day14_3") ~ 14,
                          condition %in% c("37C_day7_1", "37C_day7_2", "37C_day7_3") ~ 7,
                          condition %in% c("37C_day1_1", "37C_day1_2", "37C_day1_3") ~ 1,
                          TRUE ~ 0),
         condition = case_when(condition %in% c("37C_day14_1", "37C_day14_2", "37C_day14_3") ~ "37C day 14",
                               condition %in% c("37C_day7_1", "37C_day7_2", "37C_day7_3") ~ "37C day 7",
                               condition %in% c("37C_day1_1", "37C_day1_2", "37C_day1_3") ~ "37C day 1",
                               TRUE ~ "33C"),
         condition = factor(condition, levels = c("33C", "37C day 1", "37C day 7", "37C day 14")))



plot_df <- subset(TF_activity_sub, TF_activity_sub$source != "WT1")

ggplot(plot_df, aes(x = Time, y = score))+
  geom_point(pch = 21, alpha = 0.7, size = 2.5, aes(fill = condition))+
  facet_wrap(~ source, scales = "free_y")+
  geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, color = "gray30")+
  stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), parse = TRUE)+
  scale_fill_manual(values = colpal)+
  labs(x = "Differentiation time\n(Days)", y = "Transcription factor activity")+
  theme_few() + theme(legend.position = "none")+
  scale_x_continuous(breaks = c(0, 1, 7, 14))



plot_df2 <- subset(TF_activity_sub, TF_activity_sub$source == "WT1")
wt1_aov <- lm(score ~ condition, plot_df2)
summary(wt1_aov)

DunnettTest(score ~ condition, plot_df2, control = "33C")



ggplot(plot_df2, aes(x = condition, y = score, fill = condition))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch = 21, alpha = 0.7, dotsize = 1.5)+
  stat_summary(fun = mean, geom = "crossbar", aes(color = condition))+
  stat_compare_means(comparisons = list(c("33C", "37C day 1"), c("33C", "37C day 7"), c("33C", "37C day 14")),
                     method = "t.test", label = "p.signif")+
  scale_fill_manual(values = colpal)+
  scale_color_manual(values = colpal)+
  theme_few()+ theme(legend.position = "none")+
  labs(x =  NULL, y = "WT1 activity")



TF_activity_sub2 <- TF_activity %>%
  filter(source %in% c("E2F1", "E2F2", "E2F3", "FOXM1", "MYBL2")) %>% #Cell cycle regulators
  mutate(Time = case_when(condition %in% c("37C_day14_1", "37C_day14_2", "37C_day14_3") ~ 14,
                          condition %in% c("37C_day7_1", "37C_day7_2", "37C_day7_3") ~ 7,
                          condition %in% c("37C_day1_1", "37C_day1_2", "37C_day1_3") ~ 1,
                          TRUE ~ 0),
         condition = case_when(condition %in% c("37C_day14_1", "37C_day14_2", "37C_day14_3") ~ "37C day 14",
                               condition %in% c("37C_day7_1", "37C_day7_2", "37C_day7_3") ~ "37C day 7",
                               condition %in% c("37C_day1_1", "37C_day1_2", "37C_day1_3") ~ "37C day 1",
                               TRUE ~ "33C"),
         condition = factor(condition, levels = c("33C", "37C day 1", "37C day 7", "37C day 14")))



ggplot(TF_activity_sub2, aes(x = Time, y = score))+
  geom_point(pch = 21, alpha = 0.7, size = 2.5, aes(fill = condition))+
  facet_wrap(~ source, scales = "free_y")+
  geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, color = "gray30")+
  stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")), parse = TRUE)+
  scale_fill_manual(values = colpal)+
  labs(x = "Differentiation time\n(Days)", y = "Transcription factor activity")+
  theme_few() + theme(legend.position = "none")+
  scale_x_continuous(breaks = c(0, 1, 7, 14))



#Add on: Western blot quantification heatmap
library(readxl)
WB_quant <- read_excel("Skole/UiB/MCD-model/RNA-seq/WB_quant.xlsx", 
                       sheet = "Sheet2")

WB_quant$Time <- factor(WB_quant$Time, levels = c("33C", "37C_d1", "37C_d7", "37C_d14"))
WB_quant$Lane <- factor(WB_quant$Lane, levels = c("NPHS2", "NPHS1", "pNPHS1"))

ggplot(WB_quant, aes(x = Time, y = Lane, fill = Normalized_Int))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "#AA3377", midpoint = 0.5, name = "Normalized band intensity")+
  theme_few()+ theme(text = element_blank())

