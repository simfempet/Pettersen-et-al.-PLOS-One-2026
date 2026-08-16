library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggpointdensity)
library(RColorBrewer)
library(sf)



# Excel sheets and shapefiles
Sheets <- c("Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4",
            "PP2 trt1", "PP2 trt2", "PP2 trt3", "PP2 trt4",
            "LPS trt1", "LPS trt2", "LPS trt3", "LPS trt4",
            "Undiff 1", "Undiff 2", "Undiff 3", "Undiff 4")



file_path <- "C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Imaging/P_NPSH1.xlsx"



shp_dir <- "C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Imaging/P_NPSH1/"



#Mapping Excel sheets -> shapefiles
Shapefiles <- c(
  "Ctrl1"    = paste0(shp_dir, "Ctrl1_pseudocells.shp"),
  "Ctrl2"    = paste0(shp_dir, "Ctrl2_pseudocells.shp"),
  "Ctrl3"    = paste0(shp_dir, "Ctrl3_pseudocells.shp"),
  "Ctrl4"    = paste0(shp_dir, "Ctrl4_pseudocells.shp"),
  
  "PP2 trt1" = paste0(shp_dir, "PP2_1_pseudocells.shp"),
  "PP2 trt2" = paste0(shp_dir, "PP2_2_pseudocells.shp"),
  "PP2 trt3" = paste0(shp_dir, "PP2_3_pseudocells.shp"),
  "PP2 trt4" = paste0(shp_dir, "PP2_4_pseudocells.shp"),
  
  "LPS trt1" = paste0(shp_dir, "LPS1_pseudocells.shp"),
  "LPS trt2" = paste0(shp_dir, "LPS2_pseudocells.shp"),
  "LPS trt3" = paste0(shp_dir, "LPS3_pseudocells.shp"),
  "LPS trt4" = paste0(shp_dir, "LPS4_pseudocells.shp"),
  
  "Undiff 1" = paste0(shp_dir, "Undiff1_pseudocells.shp"),
  "Undiff 2" = paste0(shp_dir, "Undiff2_pseudocells.shp"),
  "Undiff 3" = paste0(shp_dir, "Undiff3_pseudocells.shp"),
  "Undiff 4" = paste0(shp_dir, "Undiff4_pseudocells.shp"))



data_list <- lapply(Sheets, function(s) {
  read_xlsx(file_path, sheet = s)
})
names(data_list) <- Sheets



process_condition <- function(df, shp_path, img_width = 5734, img_height = 5734) {
  
  sf_obj <- st_read(shp_path, options = c("IGNORE_MFIELDS=YES"), quiet = TRUE)
  sf_obj$FID <- sf_obj$FID + 1  
  
  scale_x <- img_width  / max(df$X, na.rm = TRUE)
  scale_y <- img_height / max(df$Y, na.rm = TRUE)
  
  df_scaled <- df %>%
    mutate(X_scaled = X * scale_x,
           Y_scaled = Y * scale_y)
  
  pts <- st_as_sf(df_scaled, coords = c("X_scaled", "Y_scaled"), crs = NA)
  
  joined <- st_join(pts, sf_obj["FID"], join = st_within)
  
  final <- joined %>%
    rename(Cell = FID) %>%
    st_drop_geometry() %>%
    filter(!is.na(Cell))
  
  return(final)
}



#Run processing function
results <- mapply(FUN = process_condition,
                  df = data_list,
                  shp_path = Shapefiles[Sheets],
                  SIMPLIFY = FALSE)



Ctrl1_final   <- results$Ctrl1
Ctrl2_final   <- results$Ctrl2
Ctrl3_final   <- results$Ctrl3
Ctrl4_final   <- results$Ctrl4



PP2trt1_final <- results$`PP2 trt1`
PP2trt2_final <- results$`PP2 trt2`
PP2trt3_final <- results$`PP2 trt3`
PP2trt4_final <- results$`PP2 trt4`



LPStrt1_final <- results$`LPS trt1`
LPStrt2_final <- results$`LPS trt2`
LPStrt3_final <- results$`LPS trt3`
LPStrt4_final <- results$`LPS trt4`



Undiff1_final <- results$`Undiff 1`
Undiff2_final <- results$`Undiff 2`
Undiff3_final <- results$`Undiff 3`
Undiff4_final <- results$`Undiff 4`



#Define Image
Ctrl1_final$ImageRegion <- 1
Ctrl2_final$ImageRegion <- 2
Ctrl3_final$ImageRegion <- 3
Ctrl4_final$ImageRegion <- 4



PP2trt1_final$ImageRegion <- 1
PP2trt2_final$ImageRegion <- 2
PP2trt3_final$ImageRegion <- 3
PP2trt4_final$ImageRegion <- 4



LPStrt1_final$ImageRegion <- 1
LPStrt2_final$ImageRegion <- 2
LPStrt3_final$ImageRegion <- 3
LPStrt4_final$ImageRegion <- 4



Undiff1_final$ImageRegion <- 1
Undiff2_final$ImageRegion <- 2
Undiff3_final$ImageRegion <- 3
Undiff4_final$ImageRegion <- 4



#Create unique cell ID
Ctrl1_final <- Ctrl1_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
Ctrl2_final <- Ctrl2_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
Ctrl3_final <- Ctrl3_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
Ctrl4_final <- Ctrl4_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))



PP2trt1_final <- PP2trt1_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
PP2trt2_final <- PP2trt2_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
PP2trt3_final <- PP2trt3_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
PP2trt4_final <- PP2trt4_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))



LPStrt1_final <- LPStrt1_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
LPStrt2_final <- LPStrt2_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
LPStrt3_final <- LPStrt3_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
LPStrt4_final <- LPStrt4_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))



Undiff1_final <- Undiff1_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
Undiff2_final <- Undiff2_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
Undiff3_final <- Undiff3_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))
Undiff4_final <- Undiff4_final %>%
  mutate(Cell_ID = paste(ImageRegion, Cell, sep = "_"))



#Combine
Ctrl_df <- rbind(Ctrl1_final, Ctrl2_final, Ctrl3_final, Ctrl4_final)



PP2_df <- rbind(PP2trt1_final, PP2trt2_final, PP2trt3_final, PP2trt4_final)



LPS_df <- rbind(LPStrt1_final, LPStrt2_final, LPStrt3_final, LPStrt4_final)



Undiff_df <- rbind(Undiff1_final, Undiff2_final, Undiff3_final, Undiff4_final)



Ctrl_df$Condition <- "Ctrl"
PP2_df$Condition <- "PP2 trt"
LPS_df$Condition <- "LPS trt"
Undiff_df$Condition <- "Undiff"



res_df <- rbind(Ctrl_df, LPS_df, PP2_df, Undiff_df)



# Compute aggregate data
res_df_summary <- res_df %>%
  group_by(Condition, Cell_ID) %>%
  summarise(mean_IntDen = mean(RawIntDen / Area, na.rm = TRUE),
            .groups = "drop")



foci_per_cell <- res_df %>%
  group_by(Condition, Cell_ID) %>%
  summarise(Foci_per_Cell = n(),
            .groups = "drop")



#Stats
#Intensity distributions
wilcox_LPS <- wilcox.test(mean_IntDen ~ Condition,
                          data = res_df_summary %>% 
                            filter(Condition %in% c("Ctrl", "LPS trt")))

wilcox_PP2 <- wilcox.test(mean_IntDen ~ Condition,
                          data = res_df_summary %>% 
                            filter(Condition %in% c("Ctrl", "PP2 trt")))

wilcox_Undiff <- wilcox.test(mean_IntDen ~ Condition,
                             data = res_df_summary %>% 
                               filter(Condition %in% c("Ctrl", "Undiff")))



wilcox_LPS
wilcox_PP2
wilcox_Undiff



pvals <- c(wilcox_LPS$p.value, wilcox_PP2$p.value, wilcox_Undiff$p.value)
pvals_adj <- p.adjust(pvals, method = "bonferroni")



#Foci per cell
wilcox_Undiff2 <- wilcox.test(Foci_per_Cell ~ Condition,
                             data = foci_per_cell %>% 
                               filter(Condition %in% c("Ctrl", "Undiff")))



# Plot
plot_df <- res_df_summary %>%
  filter(Condition != "Undiff")



ggplot(plot_df, aes(x = Condition, y = mean_IntDen, fill = Condition))+
  geom_jitter(alpha = 0.2, pch = 21)+
  stat_summary(fun = median, geom = "crossbar")+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(5, 7, 8)]))+
  stat_compare_means(comparisons=list(c("Ctrl", "PP2 trt"), c("Ctrl", "LPS trt")), method="wilcox.test", label="p.signif", 
                     tip.length=0.02, bracket.size = 0.65, size = 4.2)+
  labs(x = NULL, y = expression(atop("NPHS1-pY1176/1193", 
                                     "(Foci intensity per cell)")))+ 
  theme_few()+ theme(legend.position = "none", 
        axis.text.x = element_text(angle = 25, hjust = 1), 
        axis.title.y = element_text(hjust = 0.5))



plot_df2 <- res_df_summary %>%
  filter(Condition %in% c("Ctrl", "Undiff")) %>%
  mutate(Condition = case_when(Condition == "Ctrl" ~ "37C day 14",
                               TRUE ~ "33C"))



p1 <- ggplot(plot_df2, aes(x = Condition, y = mean_IntDen, fill = Condition))+
  geom_jitter(alpha = 0.2, pch = 21)+
  stat_summary(fun = median, geom = "crossbar")+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(1, 4)]))+
  stat_compare_means(comparisons=list(c("33C", "37C day 14")), method="wilcox.test", label="p.signif")+
  labs(x = NULL, y = "Intensity per cell", subtitle = "NPHS1-pY1176/1193")+
  ylim(0, 95000)+
  theme_few()+ theme(legend.position = "none", 
                     axis.text.x = element_text(angle = 25, hjust = 1), 
                     axis.title.y = element_text(hjust = 0.5))


plot_df3 <- foci_per_cell %>%
  filter(Condition %in% c("Ctrl", "Undiff")) %>%
  mutate(Condition = case_when(Condition == "Ctrl" ~ "37C day 14",
                               TRUE ~ "33C"))



p2 <- ggplot(plot_df3, aes(x = Condition, y = Foci_per_Cell, fill = Condition))+
  geom_jitter(alpha = 0.2, pch = 21)+
  stat_summary(fun = median, geom = "crossbar")+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(1, 4)]))+
  stat_compare_means(comparisons=list(c("33C", "37C day 14")), method="wilcox.test", label="p.signif")+
  labs(x = NULL, y = "Foci per cell", subtitle = " ")+ 
  ylim(0, 400)+
  theme_few()+ theme(legend.position = "none", 
                     axis.text.x = element_text(angle = 25, hjust = 1), 
                     axis.title.y = element_text(hjust = 0.5))



ggarrange(p1, p2, align = "h")


# Plotting shapefiles
Ctrl1_sf <- st_read(
  "C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Imaging/P_NPSH1/Ctrl1_pseudocells.shp",
  options = c("IGNORE_MFIELDS=YES"))

Ctrl1_sf$FID <- Ctrl1_sf$FID + 1

img_height <- 5734
img_width  <- 5734

scale_x <- img_width  / max(Ctrl1_final$X)
scale_y <- img_height / max(Ctrl1_final$Y)


Ctrl1_final <- Ctrl1_final %>%
  mutate(X_scaled = X * scale_x, Y_scaled = Y * scale_y)

Ctrl1_pts <- st_as_sf(Ctrl1_final, coords = c("X_scaled", "Y_scaled"), crs = NA)


a1 <- ggplot() +
  geom_sf(data = Ctrl1_sf, aes(fill = as.numeric(FID)), size = 0.3) +
  geom_sf(data = Ctrl1_pts, color = "cyan", size = 1.2, alpha = 0.3) +
  scale_fill_viridis_c(option = "rocket")+ xlim(2000, 4500)+ ylim(2000, 4500)+
  theme_void()+ theme(legend.position = "none")+ labs(subtitle = "37C day 14")



Undiff1_sf <- st_read(
  "C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Imaging/P_NPSH1/Undiff1_pseudocells.shp",
  options = c("IGNORE_MFIELDS=YES"))

Undiff1_sf$FID <- Undiff1_sf$FID + 1

img_height <- 5734
img_width  <- 5734

scale_x <- img_width  / max(Undiff1_final$X)
scale_y <- img_height / max(Undiff1_final$Y)


Undiff1_final <- Undiff1_final %>%
  mutate(X_scaled = X * scale_x, Y_scaled = Y * scale_y)

Undiff1_pts <- st_as_sf(Undiff1_final, coords = c("X_scaled", "Y_scaled"), crs = NA)


a2 <- ggplot() +
  geom_sf(data = Undiff1_sf, aes(fill = as.numeric(FID)), size = 0.3) +
  geom_sf(data = Undiff1_pts, color = "cyan", size = 1.2, alpha = 0.3) +
  scale_fill_viridis_c(option = "rocket")+ xlim(2000, 4500)+ ylim(2000, 4500)+
  theme_void()+ theme(legend.position = "none")+ labs(subtitle = "33C")



ggarrange(a2, a1, align = "h")



#Save outputs
write.csv(res_df_summary, file = "C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Imaging/pNPHS1_intensity.csv")
write.csv(foci_per_cell, file = "C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/Imaging/pNPHS1_counts.csv")