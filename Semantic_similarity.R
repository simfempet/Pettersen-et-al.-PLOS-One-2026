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



#In vitro results
GO_BP_LPS_up <- readRDS("C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/RNA-seq/GO_LPS_up.RDS")
GO_BP_LPS_down <- readRDS("C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/RNA-seq/GO_LPS_down.RDS")



GO_BP_PP2_up <- readRDS("C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/RNA-seq/GO_PP2_up.RDS")
GO_BP_PP2_down <- readRDS("C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/RNA-seq/GO_PP2_down.RDS")



#Helper function for weighted similarity
combine_up_down_similarity <- function(up_exp1_ids, up_exp2_ids,
                                       down_exp1_ids, down_exp2_ids,
                                       S_up, S_down) {
  n_up   <- length(unique(c(up_exp1_ids, up_exp2_ids)))
  n_down <- length(unique(c(down_exp1_ids, down_exp2_ids)))
  
  if (is.na(S_up) && is.na(S_down)) {
    return(NA_real_)
  }
  if (is.na(S_up)) return(S_down)
  if (is.na(S_down)) return(S_up)
  
  overall_weighted <- (n_up * S_up + n_down * S_down) / (n_up + n_down)
  
  overall_mean <- mean(c(S_up, S_down))
  
  list(
    n_up = n_up, n_down = n_down,
    S_up = S_up, S_down = S_down,
    overall_weighted = overall_weighted,
    overall_unweighted = overall_mean
  )
}




#Diabetic nephropathy (GSE142025)
dat <- GSE142025_raw_counts_GRCh38.p13_NCBI
dat <- dat[!duplicated(dat$GeneID), ]
dat <- dat[!(dat$GeneID == ""), ]
dat <- dat[!is.na(dat$GeneID), ]



rownames(dat) <- dat$GeneID
dat <- dat[, -1]



gse <- getGEO(GEO = "GSE142025", GSEMatrix = T)



metadata <- pData(phenoData(gse[[1]]))



colData <- metadata[,c(1, 2, 37)]
colnames(colData)[3] <- "Condition"



all(colnames(dat) %in% rownames(colData))
all(colnames(dat) == rownames(colData))



dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~ Condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds$Condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)



res <- results(dds, contrast = c("Condition", "Advanced_DN", "Control"), alpha = 0.05)
res_df <- data.frame(res)
res_df <- res_df %>% 
  mutate(Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up",
                                log2FoldChange < 0 & padj <= 0.05 ~ "Down",
                                TRUE ~ "Unchanged"),
         GeneID = rownames(res_df))
res_df <- res_df[! is.na(res_df$Expression), ]



ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id"),
  filters = "entrezgene_id",
  values = res_df$GeneID,
  mart = ensembl)



res_df <- merge(res_df, gene_map, by.x = "GeneID", by.y = "entrezgene_id", all.x = TRUE)
res_df <- res_df[res_df$hgnc_symbol != "", ]
res_df <- res_df[! is.na(res_df$Expression), ]



DN_up <- (subset(res_df, res_df$Expression == "Up"))$hgnc_symbol
DN_down <- (subset(res_df, res_df$Expression == "Down"))$hgnc_symbol



GO_BP_DN_up <- enrichGO(gene = DN_up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$hgnc_symbol)
GO_BP_DN_up_df <- as.data.frame(GO_BP_DN_up)



GO_BP_DN_down <- enrichGO(gene = DN_down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$hgnc_symbol)
GO_BP_DN_down_df <- as.data.frame(GO_BP_DN_down)



hsGO <- godata('org.Hs.eg.db', ont="BP")



sim_matrix_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_DN_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_up[is.na(sim_matrix_LPS_up)] <- 0
overall_sim_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_DN_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.778
LPSvDN_up <- overall_sim_LPS_up



sim_matrix_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_DN_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_down[is.na(sim_matrix_LPS_down)] <- 0
overall_sim_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_DN_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.439
LPSvDN_down <- overall_sim_LPS_down



sim_matrix_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_DN_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_up[is.na(sim_matrix_PP2_up)] <- 0
overall_sim_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_DN_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.317
PP2vDN_up <- overall_sim_PP2_up



sim_matrix_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_DN_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_down[is.na(sim_matrix_PP2_down)] <- 0
overall_sim_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_DN_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.462
PP2vDN_down <- overall_sim_PP2_down



LPS_combined_similarity_DN <- combine_up_down_similarity(up_exp1_ids   = GO_BP_LPS_up$ID,
                                                      up_exp2_ids   = GO_BP_DN_up$ID,
                                                      down_exp1_ids = GO_BP_LPS_down$ID,
                                                      down_exp2_ids = GO_BP_DN_down$ID,
                                                      S_up = overall_sim_LPS_up,
                                                      S_down = overall_sim_LPS_down)
LPS_combined_similarity_DN$overall_weighted #0.667
LPS_combined_similarity_DN$overall_unweighted #0.608



PP2_combined_similarity_DN <- combine_up_down_similarity(up_exp1_ids   = GO_BP_PP2_up$ID,
                                                      up_exp2_ids   = GO_BP_DN_up$ID,
                                                      down_exp1_ids = GO_BP_PP2_down$ID,
                                                      down_exp2_ids = GO_BP_DN_down$ID,
                                                      S_up = overall_sim_PP2_up,
                                                      S_down = overall_sim_PP2_down)
PP2_combined_similarity_DN$overall_weighted #0.378
PP2_combined_similarity_DN$overall_unweighted # 0.389



colpal <- c(brewer.pal(12, "Paired")[7], "hotpink1")
thr_lower <- 0.9
thr_upper <- 0.999



LPS_ids <- GO_BP_LPS_up$ID
DN_ids  <- GO_BP_DN_up$ID
LPS_desc <- GO_BP_LPS_up$Description
DN_desc  <- GO_BP_DN_up$Description



vertices <- data.frame(
  name = c(paste0(LPS_ids, "_LPS"), paste0(DN_ids, "_DN")),
  label = c(LPS_desc, DN_desc),      
  source = c(rep("LPS_up", length(LPS_ids)),
             rep("DN_up", length(DN_ids))),
  stringsAsFactors = FALSE
)



edges <- expand.grid(from = LPS_ids, to = DN_ids, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(weight = sim_matrix_LPS_up[from, to]) %>%
  ungroup() %>%
  filter(weight >= thr_lower & weight < thr_upper) %>%
  mutate(
    from = paste0(from, "_LPS"),
    to   = paste0(to, "_DN")
  )


g_DN <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)
g_DN <- delete_vertices(g_DN, V(g_DN)[degree(g_DN) == 0])
V(g_DN)$desc <- V(g_DN)$label



pattern <- "immune|inflammatory|cytokine|chemokine"
V(g_DN)$label <- ""
immune_nodes <- str_detect(V(g_DN)$desc, regex(pattern, ignore_case = TRUE))
V(g_DN)$label[immune_nodes] <- V(g_DN)$desc[immune_nodes]



tg_DN <- as_tbl_graph(g_DN)
a1 <- ggraph(tg_DN, layout = "linear", circular = T) +
  geom_edge_link(aes(width = weight), color = "gray55", alpha = 0.2) +
  geom_node_point(aes(color = source), size = 3, alpha = 0.7) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold",
                 color = "black", bg.color = "white", bg.r = 0.15) +
  scale_color_manual(values = c("LPS_up" = colpal[1], "DN_up" = colpal[2])) +
  theme_void()



#Hypertensive nephropathy (GSE166239)
dat <- GSE166239_Nordbo_et_al_counts
dat <- dat[!duplicated(dat$contig), ]
dat <- dat[!(dat$contig == ""), ]
dat <- dat[!is.na(dat$contig), ]



rownames(dat) <- dat$contig
dat <- dat[, -1]



gse <- getGEO(GEO = "GSE166239", GSEMatrix = T)



metadata <- pData(phenoData(gse[[1]]))
colData <- metadata[,c(1, 2)]



colData$Condition <- c(rep("Control", 6), rep("HN", 6), rep("DN", 6))
rownames(colData) <-  colData$title
dat <- dat[, rownames(colData)]



all(colnames(dat) %in% rownames(colData))
all(colnames(dat) == rownames(colData))



dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~ Condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds$Condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)



res <- results(dds, contrast = c("Condition", "HN", "Control"), alpha = 0.05)
res_df <- data.frame(res)
res_df <- res_df %>% 
  mutate(Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up",
                                log2FoldChange < 0 & padj <= 0.05 ~ "Down",
                                TRUE ~ "Unchanged"),
         GeneID = rownames(res_df))
res_df <- res_df[! is.na(res_df$Expression), ]



gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res_df$GeneID,
  mart = ensembl)



res_df <- merge(res_df, gene_map, by.x = "GeneID", by.y = "ensembl_gene_id", all.x = TRUE)
res_df <- res_df[res_df$hgnc_symbol != "", ]
res_df <- res_df[! is.na(res_df$Expression), ]



HN_up <- (subset(res_df, res_df$Expression == "Up"))$hgnc_symbol
HN_down <- (subset(res_df, res_df$Expression == "Down"))$hgnc_symbol



GO_BP_HN_up <- enrichGO(gene = HN_up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$hgnc_symbol)
GO_BP_HN_up_df <- as.data.frame(GO_BP_HN_up)



GO_BP_HN_down <- enrichGO(gene = HN_down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$hgnc_symbol)
GO_BP_HN_down_df <- as.data.frame(GO_BP_HN_down)



sim_matrix_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_HN_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_up[is.na(sim_matrix_LPS_up)] <- 0
overall_sim_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_HN_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.787
LPSvHN_up <- overall_sim_LPS_up



sim_matrix_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_HN_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_down[is.na(sim_matrix_LPS_down)] <- 0
overall_sim_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_HN_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.346
LPSvHN_down <- overall_sim_LPS_down



sim_matrix_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_HN_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_up[is.na(sim_matrix_PP2_up)] <- 0
overall_sim_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_HN_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.331
PP2vHN_up <- overall_sim_PP2_up



sim_matrix_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_HN_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_down[is.na(sim_matrix_PP2_down)] <- 0
overall_sim_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_HN_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.425
PP2vHN_down <- overall_sim_PP2_down



LPS_combined_similarity_HN <- combine_up_down_similarity(up_exp1_ids   = GO_BP_LPS_up$ID,
                                                      up_exp2_ids   = GO_BP_HN_up$ID,
                                                      down_exp1_ids = GO_BP_LPS_down$ID,
                                                      down_exp2_ids = GO_BP_HN_down$ID,
                                                      S_up = overall_sim_LPS_up,
                                                      S_down = overall_sim_LPS_down)
LPS_combined_similarity_HN$overall_weighted#0.691
LPS_combined_similarity_HN$overall_unweighted#0.566



PP2_combined_similarity_HN <- combine_up_down_similarity(up_exp1_ids   = GO_BP_PP2_up$ID,
                                                      up_exp2_ids   = GO_BP_HN_up$ID,
                                                      down_exp1_ids = GO_BP_PP2_down$ID,
                                                      down_exp2_ids = GO_BP_HN_down$ID,
                                                      S_up = overall_sim_PP2_up,
                                                      S_down = overall_sim_PP2_down)
PP2_combined_similarity_HN$overall_weighted#0.356
PP2_combined_similarity_HN$overall_unweighted#0.378



colpal <- c(brewer.pal(12, "Paired")[7], "darkseagreen2")



LPS_ids <- GO_BP_LPS_up$ID
HN_ids  <- GO_BP_HN_up$ID
LPS_desc <- GO_BP_LPS_up$Description
HN_desc  <- GO_BP_HN_up$Description



vertices <- data.frame(
  name = c(paste0(LPS_ids, "_LPS"), paste0(HN_ids, "_HN")),
  label = c(LPS_desc, HN_desc),      
  source = c(rep("LPS_up", length(LPS_ids)),
             rep("HN_up", length(HN_ids))),
  stringsAsFactors = FALSE
)



edges <- expand.grid(from = LPS_ids, to = HN_ids, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(weight = sim_matrix_LPS_up[from, to]) %>%
  ungroup() %>%
  filter(weight >= thr_lower & weight < thr_upper) %>%
  mutate(
    from = paste0(from, "_LPS"),
    to   = paste0(to, "_HN")
  )


g_HN <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)
g_HN <- delete_vertices(g_HN, V(g_HN)[degree(g_HN) == 0])
V(g_HN)$desc <- V(g_HN)$label



V(g_HN)$label <- ""
immune_nodes <- str_detect(V(g_HN)$desc, regex(pattern, ignore_case = TRUE))
V(g_HN)$label[immune_nodes] <- V(g_HN)$desc[immune_nodes]



tg_HN <- as_tbl_graph(g_HN)
a2 <- ggraph(tg_HN, layout = "linear", circular = T) +
  geom_edge_link(aes(width = weight), color = "gray55", alpha = 0.2) +
  geom_node_point(aes(color = source), size = 3, alpha = 0.7) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold",
                 color = "black", bg.color = "white", bg.r = 0.15) +
  scale_color_manual(values = c("LPS_up" = colpal[1], "HN_up" = colpal[2])) +
  theme_void()



#IgA nephropathy (GSE141295)
dat <- GSE141295_raw_counts_GRCh38.p13_NCBI
dat <- dat[!duplicated(dat$GeneID), ]
dat <- dat[!(dat$GeneID == ""), ]
dat <- dat[!is.na(dat$GeneID), ]



rownames(dat) <- dat$GeneID
dat <- dat[, -1]



gse <- getGEO(GEO = "GSE141295", GSEMatrix = T)



metadata <- pData(phenoData(gse[[1]]))
colData <- metadata[, c(2, 40)]



colData$Condition <- c(rep("IgAN", 15), rep("Control", 9))



all(colnames(dat) %in% rownames(colData))
all(colnames(dat) == rownames(colData))



dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData, design = ~ Condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds$Condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)



res <- results(dds, contrast = c("Condition", "IgAN", "Control"), alpha = 0.05)
res_df <- data.frame(res)
res_df <- res_df %>% 
  mutate(Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up",
                                log2FoldChange < 0 & padj <= 0.05 ~ "Down",
                                TRUE ~ "Unchanged"),
         GeneID = rownames(res_df))



gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id"),
  filters = "entrezgene_id",
  values = res_df$GeneID,
  mart = ensembl)



res_df <- merge(res_df, gene_map, by.x = "GeneID", by.y = "entrezgene_id", all.x = TRUE)
res_df <- res_df[res_df$hgnc_symbol != "", ]
res_df <- res_df[! is.na(res_df$Expression), ]



IgAN_up <- (subset(res_df, res_df$Expression == "Up"))$hgnc_symbol
IgAN_down <- (subset(res_df, res_df$Expression == "Down"))$hgnc_symbol



GO_BP_IgAN_up <- enrichGO(gene = IgAN_up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$hgnc_symbol)
GO_BP_IgAN_up_df <- as.data.frame(GO_BP_IgAN_up)



GO_BP_IgAN_down <- enrichGO(gene = IgAN_down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$hgnc_symbol)
GO_BP_IgAN_down_df <- as.data.frame(GO_BP_IgAN_down)



sim_matrix_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_IgAN_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_up[is.na(sim_matrix_LPS_up)] <- 0
overall_sim_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_IgAN_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.744
LPSvIgAN_up <- overall_sim_LPS_up



sim_matrix_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_IgAN_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_down[is.na(sim_matrix_LPS_down)] <- 0
overall_sim_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_IgAN_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.197
LPSvIgAN_down <- overall_sim_LPS_down



sim_matrix_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_IgAN_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_up[is.na(sim_matrix_PP2_up)] <- 0
overall_sim_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_IgAN_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.348
PP2vIgAN_up <- overall_sim_PP2_up



sim_matrix_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_IgAN_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_down[is.na(sim_matrix_PP2_down)] <- 0
overall_sim_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_IgAN_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.26
PP2vIgAN_down <- overall_sim_PP2_down



LPS_combined_similarity_IgAN <- combine_up_down_similarity(up_exp1_ids   = GO_BP_LPS_up$ID,
                                                      up_exp2_ids   = GO_BP_IgAN_up$ID,
                                                      down_exp1_ids = GO_BP_LPS_down$ID,
                                                      down_exp2_ids = GO_BP_IgAN_down$ID,
                                                      S_up = overall_sim_LPS_up,
                                                      S_down = overall_sim_LPS_down)
LPS_combined_similarity_IgAN$overall_weighted#0.686
LPS_combined_similarity_IgAN$overall_unweighted#0.470



PP2_combined_similarity_IgAN <- combine_up_down_similarity(up_exp1_ids   = GO_BP_PP2_up$ID,
                                                      up_exp2_ids   = GO_BP_IgAN_up$ID,
                                                      down_exp1_ids = GO_BP_PP2_down$ID,
                                                      down_exp2_ids = GO_BP_IgAN_down$ID,
                                                      S_up = overall_sim_PP2_up,
                                                      S_down = overall_sim_PP2_down)
PP2_combined_similarity_IgAN$overall_weighted#0.335
PP2_combined_similarity_IgAN$overall_unweighted#0.304



colpal <- c(brewer.pal(12, "Paired")[7], "steelblue")



LPS_ids <- GO_BP_LPS_up$ID
IgAN_ids  <- GO_BP_IgAN_up$ID
LPS_desc <- GO_BP_LPS_up$Description
IgAN_desc  <- GO_BP_IgAN_up$Description



vertices <- data.frame(
  name = c(paste0(LPS_ids, "_LPS"), paste0(IgAN_ids, "_IgAN")),
  label = c(LPS_desc, IgAN_desc),      
  source = c(rep("LPS_up", length(LPS_ids)),
             rep("IgAN_up", length(IgAN_ids))),
  stringsAsFactors = FALSE
)



edges <- expand.grid(from = LPS_ids, to = IgAN_ids, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(weight = sim_matrix_LPS_up[from, to]) %>%
  ungroup() %>%
  filter(weight >= thr_lower & weight < thr_upper) %>%
  mutate(
    from = paste0(from, "_LPS"),
    to   = paste0(to, "_IgAN")
  )



g_IgAN <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)
g_IgAN <- delete_vertices(g_IgAN, V(g_IgAN)[degree(g_IgAN) == 0])
V(g_IgAN)$desc <- V(g_IgAN)$label



V(g_IgAN)$label <- ""
immune_nodes <- str_detect(V(g_IgAN)$desc, regex(pattern, ignore_case = TRUE))
V(g_IgAN)$label[immune_nodes] <- V(g_IgAN)$desc[immune_nodes]



tg_IgAN <- as_tbl_graph(g_IgAN)
a3 <- ggraph(tg_IgAN, layout = "linear", circular = T) +
  geom_edge_link(aes(width = weight), color = "gray55", alpha = 0.2) +
  geom_node_point(aes(color = source), size = 3, alpha = 0.7) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold",
                 color = "black", bg.color = "white", bg.r = 0.15) +
  scale_color_manual(values = c("LPS_up" = colpal[1], "IgAN_up" = colpal[2])) +
  theme_void()



#MCD (GSE216841)
dat <- GSE216841_MNd_ncounts_annot
dat <- dat[!duplicated(dat$hgnc_GRCh38p12), ]
dat <- dat[!(dat$hgnc_GRCh38p12 == ""), ]
dat <- dat[!is.na(dat$hgnc_GRCh38p12), ]



counts <- dat[, c(5:38)]
counts <- round(counts, digits = 0)
rownames(counts) <- dat$hgnc_GRCh38p12



gse <- getGEO(GEO = "GSE216841", GSEMatrix = T)



metadata <- pData(phenoData(gse[[1]]))
colData <- metadata[,c(1, 11)]
rownames(colData) <- colnames(counts)
colData$Type <- c(rep("NC", 8), rep("MN", 12), rep("MCD", 14))



all(colnames(counts) %in% rownames(colData))
all(colnames(counts) == rownames(colData))



dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Type)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds$Type <- relevel(dds$Type, ref = "NC")
dds <- DESeq(dds)



res <- results(dds, contrast = c("Type", "MCD", "NC"))
res_df <- data.frame(res)
res_df <- res_df %>% 
  mutate(Expression = case_when(log2FoldChange > 0 & padj <= 0.05 ~ "Up",
                                log2FoldChange < 0 & padj <= 0.05 ~ "Down",
                                TRUE ~ "Unchanged"),
         GeneID = rownames(res_df))
res_df <- res_df[! is.na(res_df$Expression), ]



MCD_up <- (subset(res_df, res_df$Expression == "Up"))$GeneID
MCD_down <- (subset(res_df, res_df$Expression == "Down"))$GeneID



GO_BP_MCD_up <- enrichGO(gene = MCD_up, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$GeneID)
GO_BP_MCD_up_df <- as.data.frame(GO_BP_MCD_up)



GO_BP_MCD_down <- enrichGO(gene = MCD_down, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", universe = res_df$GeneID)
GO_BP_MCD_down_df <- as.data.frame(GO_BP_MCD_down)



sim_matrix_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_MCD_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_up[is.na(sim_matrix_LPS_up)] <- 0
overall_sim_LPS_up <- mgoSim(GO_BP_LPS_up$ID, GO_BP_MCD_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.751
LPSvMCD_up <- overall_sim_LPS_up



sim_matrix_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_MCD_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_LPS_down[is.na(sim_matrix_LPS_down)] <- 0
overall_sim_LPS_down <- mgoSim(GO_BP_LPS_down$ID, GO_BP_MCD_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.333
LPSvMCD_down <- overall_sim_LPS_down



sim_matrix_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_MCD_up$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_up[is.na(sim_matrix_PP2_up)] <- 0
overall_sim_PP2_up <- mgoSim(GO_BP_PP2_up$ID, GO_BP_MCD_up$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.349
PP2vMCD_up <- overall_sim_PP2_up



sim_matrix_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_MCD_down$ID, semData = hsGO, measure = "Wang", combine = NULL)
sim_matrix_PP2_down[is.na(sim_matrix_PP2_down)] <- 0
overall_sim_PP2_down <- mgoSim(GO_BP_PP2_down$ID, GO_BP_MCD_down$ID, semData = hsGO, measure = "Wang", combine = "BMA")#0.377
PP2vMCD_down <- overall_sim_PP2_down



LPS_combined_similarity_MCD <- combine_up_down_similarity(up_exp1_ids   = GO_BP_LPS_up$ID,
                                                      up_exp2_ids   = GO_BP_MCD_up$ID,
                                                      down_exp1_ids = GO_BP_LPS_down$ID,
                                                      down_exp2_ids = GO_BP_MCD_down$ID,
                                                      S_up = overall_sim_LPS_up,
                                                      S_down = overall_sim_LPS_down)
LPS_combined_similarity_MCD$overall_weighted#0.661
LPS_combined_similarity_MCD$overall_unweighted#0.542



PP2_combined_similarity_MCD <- combine_up_down_similarity(up_exp1_ids   = GO_BP_PP2_up$ID,
                                                      up_exp2_ids   = GO_BP_MCD_up$ID,
                                                      down_exp1_ids = GO_BP_PP2_down$ID,
                                                      down_exp2_ids = GO_BP_MCD_down$ID,
                                                      S_up = overall_sim_PP2_up,
                                                      S_down = overall_sim_PP2_down)
PP2_combined_similarity_MCD$overall_weighted#0.357
PP2_combined_similarity_MCD$overall_unweighted#0.363



colpal <- c(brewer.pal(12, "Paired")[7], "gray")



LPS_ids <- GO_BP_LPS_up$ID
MCD_ids  <- GO_BP_MCD_up$ID
LPS_desc <- GO_BP_LPS_up$Description
MCD_desc  <- GO_BP_MCD_up$Description



vertices <- data.frame(
  name = c(paste0(LPS_ids, "_LPS"), paste0(MCD_ids, "_MCD")),
  label = c(LPS_desc, MCD_desc),      
  source = c(rep("LPS_up", length(LPS_ids)),
             rep("MCD_up", length(MCD_ids))),
  stringsAsFactors = FALSE
)



edges <- expand.grid(from = LPS_ids, to = MCD_ids, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(weight = sim_matrix_LPS_up[from, to]) %>%
  ungroup() %>%
  filter(weight >= thr_lower & weight < thr_upper) %>%
  mutate(
    from = paste0(from, "_LPS"),
    to   = paste0(to, "_MCD")
  )



g_MCD <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)
g_MCD <- delete_vertices(g_MCD, V(g_MCD)[degree(g_MCD) == 0])
V(g_MCD)$desc <- V(g_MCD)$label



V(g_MCD)$label <- ""
immune_nodes <- str_detect(V(g_MCD)$desc, regex(pattern, ignore_case = TRUE))
V(g_MCD)$label[immune_nodes] <- V(g_MCD)$desc[immune_nodes]



tg_MCD <- as_tbl_graph(g_MCD)
a4 <- ggraph(tg_MCD, layout = "linear", circular = T) +
  geom_edge_link(aes(width = weight), color = "gray55", alpha = 0.2) +
  geom_node_point(aes(color = source), size = 3, alpha = 0.7) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold",
                 color = "black", bg.color = "white", bg.r = 0.15) +
  scale_color_manual(values = c("LPS_up" = colpal[1], "MCD_up" = colpal[2])) +
  theme_void()



#Plots
a1 <- a1+labs(subtitle = "DN")
a2 <- a2+labs(subtitle = "HN")
a3 <- a3+labs(subtitle = "IgAN")
a4 <- a4+labs(subtitle = "MCD")



ggarrange(a1, a2, a3, a4, ncol = 2, nrow = 2, align = "hv")



DN_LPS_up <- V(g_DN)$label
DN_LPS_up <- DN_LPS_up[DN_LPS_up != ""]



HN_LPS_up <- V(g_HN)$label
HN_LPS_up <- HN_LPS_up[HN_LPS_up != ""]



IgAN_LPS_up <- V(g_IgAN)$label
IgAN_LPS_up <- IgAN_LPS_up[IgAN_LPS_up != ""]



MCD_LPS_up <- V(g_MCD)$label
MCD_LPS_up <- MCD_LPS_up[MCD_LPS_up != ""]



HighSim <- list(DN_LPS_up = DN_LPS_up,
                HN_LPS_up = HN_LPS_up,
                IgAN_LPS_up = IgAN_LPS_up,
                MCD_LPS_up = MCD_LPS_up)



m1 = make_comb_mat(HighSim)
p1 <- UpSet(m1, column_title = "Immune related processes\n(semantic similarity > 0.9)",
            top_annotation =  upset_top_annotation(m1, add_numbers = T,
                                                   gp = gpar(fill = "gray95", color = "black")),
            right_annotation = upset_right_annotation(m1, add_numbers = T,
                                                      gp = gpar(fill = "gray95", color = "black")))
p1 <- grid.grabExpr(draw(p1))
p1 <- as_ggplot(p1)



plot_df <- data.frame(Treatment = rep(c("LPS", "PP2"),16),
                        Disease = rep(c(rep("DN", 2), rep("HN", 2), rep("IgAN", 2), rep("MCD", 2))),
                        Direction = c(rep("Upregulated", 8), rep("Downregulated", 8)),
                        SemSim = c(LPSvDN_up, PP2vDN_up,
                                   LPSvHN_up, PP2vHN_up,
                                   LPSvIgAN_up, PP2vIgAN_up,
                                   LPSvMCD_up, PP2vMCD_up,
                                   LPSvDN_down, PP2vDN_down,
                                   LPSvHN_down, PP2vHN_down,
                                   LPSvIgAN_down, PP2vIgAN_down,
                                   LPSvMCD_down, PP2vMCD_down))



plot_df$Direction <- factor(plot_df$Direction, levels = c("Upregulated", "Downregulated"))



p1 <- ggplot(plot_df, aes(y = Disease, x = SemSim, group = Treatment, color = Treatment, fill = Treatment))+
  geom_linerange(aes(xmin = 0, xmax = SemSim), position = position_dodge(0.5), linewidth = 1)+
  geom_point(size = 4, position = position_dodge(0.5), pch = 21, color = "black")+
  facet_wrap(~ Direction, nrow = 2)+
  labs(y = "Disease category", x = "Semantic similarity\n(biological process)")+
  scale_color_manual(values = (brewer.pal(12, "Paired")[c(7,8)]))+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(7,8)]))+
  theme_few()+ xlim(0,1)+ geom_vline(xintercept = 0.6, linetype = "dashed")



plot_df2 <- data.frame(Treatment = rep(c("LPS", "PP2"),8),
                        Disease = c(rep("DN", 2), rep("HN", 2), rep("IgAN", 2), rep("MCD", 2)),
                        SemSim = c(LPS_combined_similarity_DN$overall_weighted, PP2_combined_similarity_DN$overall_weighted,
                                   LPS_combined_similarity_HN$overall_weighted, PP2_combined_similarity_HN$overall_weighted,
                                   LPS_combined_similarity_IgAN$overall_weighted, PP2_combined_similarity_IgAN$overall_weighted,
                                   LPS_combined_similarity_MCD$overall_weighted, PP2_combined_similarity_MCD$overall_weighted))



p2 <- ggplot(plot_df2, aes(y = Disease, x = SemSim, group = Treatment, fill = Treatment, color = Treatment))+
  geom_linerange(aes(xmin = 0, xmax = SemSim), position = position_dodge(0.5), linewidth = 1)+
  geom_point(size = 4, position = position_dodge(0.5), pch = 21, color = "black")+
  labs(y = "Disease category", x = "Weighted mean\n(upregulated & downregulated)")+
  scale_color_manual(values = (brewer.pal(12, "Paired")[c(7,8)]))+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(7,8)]))+
  theme_few()+ xlim(0,1)+ geom_vline(xintercept = 0.6, linetype = "dashed")



plot_df3 <- data.frame(Values = seq(0, 1, 0.01), Rank = seq(0, 1, 0.01))



p3 <- ggplot(plot_df3, aes(x = Rank, y = "All", fill = Values, colour = Values))+
  geom_tile()+
  scale_fill_viridis_c(option = "rocket", breaks = c(0, 0.2, 0.4, 0.6, 0.75, 0.9, 1), name = "Semantic\nsimilarity")+
  scale_color_viridis_c(option = "rocket", breaks = c(0, 0.2, 0.4, 0.6, 0.75, 0.9, 1), name = "Semantic\nsimilarity")+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.75, 0.9, 1), labels = c("None (0)", "Low (0.2)", "Moderate (0.4)", "High (0.6)", "Very high (0.75)", "Near identical (0.9)", "Identical (1)"))+
  labs(x = NULL, y = NULL, subtitle = "Semantic\nsimilarity")+
  theme_few()+ theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")



ggarrange(p3, p1, p2, nrow = 3, heights = c(0.3, 1, 0.5), align = "v")


