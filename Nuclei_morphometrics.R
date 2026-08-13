library(readxl)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(ggthemes)



#LPS
Nuclei_data_LPS <- read_excel("Skole/UiB/MCD-model/RNA-seq/Nuclei_data.xlsx",
                          sheet = "Sheet1")



Nuclei_data_LPS$Cond <- factor(Nuclei_data_LPS$Cond)



Nuclei_data_LPS <- Nuclei_data_LPS %>%
  mutate(IntCorr = RawIntDen/Area) %>%
  mutate(Cond = case_when(Cond == "Ctrl_LPS" ~ "Ctrl LPS",
                          TRUE ~ "LPS (40 µg/mL)"))



LPS_df <- Nuclei_data_LPS %>%
  group_by(Cond, Sample) %>%
  summarise(n_nuclei = n(),
            across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")



LPS_lm_circ <- lm(Circ ~ Cond, LPS_df)
summary(LPS_lm_circ)



LPS_lm_sol <- lm(Solidity ~ Cond, LPS_df)
summary(LPS_lm_sol)



LPS_lm_AR <- lm(AR ~ Cond, LPS_df)
summary(LPS_lm_AR)



LPS_pvals <- c(summary(LPS_lm_circ)[[4]][2, 4],
               summary(LPS_lm_sol)[[4]][2, 4],
               summary(LPS_lm_AR)[[4]][2, 4])

LPS_pvals_adj <- p.adjust(LPS_pvals, method = "bonferroni")



#PP2
Nuclei_data_PP2 <- read_excel("Skole/UiB/MCD-model/RNA-seq/Nuclei_data.xlsx",
                              sheet = "Sheet2")



Nuclei_data_PP2$Cond <- factor(Nuclei_data_PP2$Cond)



Nuclei_data_PP2 <- Nuclei_data_PP2 %>%
  mutate(IntCorr = RawIntDen/Area) %>%
  mutate(Cond = case_when(Cond == "Ctrl_PP2'" ~ "Ctrl PP2",
                          TRUE ~ "PP2 (75 µM)"))



PP2_df <- Nuclei_data_PP2 %>%
  group_by(Cond, Sample) %>%
  summarise(n_nuclei = n(),
            across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")



PP2_lm_circ <- lm(Circ ~ Cond, PP2_df)
summary(PP2_lm_circ)



PP2_lm_sol <- lm(Solidity ~ Cond, PP2_df)
summary(PP2_lm_sol)



PP2_lm_AR <- lm(AR ~ Cond, PP2_df)
summary(PP2_lm_AR)



PP2_pvals <- c(summary(PP2_lm_circ)[[4]][2, 4],
               summary(PP2_lm_sol)[[4]][2, 4],
               summary(PP2_lm_AR)[[4]][2, 4])

PP2_pvals_adj <- p.adjust(PP2_pvals, method = "bonferroni")



#Plot
df <- rbind(LPS_df, PP2_df)
df$Cond <- factor(df$Cond, levels = c("Ctrl LPS", "LPS (40 µg/mL)",
                                      "Ctrl PP2", "PP2 (75 µM)"))

p1 <- ggplot(df, aes(x = Cond, y = Circ, fill = Cond))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch = 21, color = "black", alpha = 0.5)+ ylim(0, 1)+
  stat_summary(fun = mean, geom = "crossbar", aes(color = Cond))+
  stat_compare_means(comparisons = list(c("Ctrl LPS", "LPS (40 µg/mL)"), c("Ctrl PP2", "PP2 (75 µM)")), label = "p.signif")+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(5,7,6,8)]))+
  scale_color_manual(values = (brewer.pal(12, "Paired")[c(5,7,6,8)]))+
  theme_few()+ theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))+
  labs(x = NULL, y = "Nuclei circularity")

p2 <- ggplot(df, aes(x = Cond, y = Solidity, fill = Cond))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch = 21, color = "black", alpha = 0.5)+ ylim(0, 1.5)+
  stat_summary(fun = mean, geom = "crossbar", aes(color = Cond))+
  stat_compare_means(comparisons = list(c("Ctrl LPS", "LPS (40 µg/mL)"), c("Ctrl PP2", "PP2 (75 µM)")), label = "p.signif")+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(5,7,6,8)]))+
  scale_color_manual(values = (brewer.pal(12, "Paired")[c(5,7,6,8)]))+
  theme_few()+ theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))+
  labs(x = NULL, y = "Nuclei solidity")

p3 <- ggplot(df, aes(x = Cond, y = AR, fill = Cond))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch = 21, color = "black", alpha = 0.5)+ylim(0, 2)+
  stat_summary(fun = mean, geom = "crossbar", aes(color = Cond))+
  stat_compare_means(comparisons = list(c("Ctrl LPS", "LPS (40 µg/mL)"), c("Ctrl PP2", "PP2 (75 µM)")), label = "p.signif")+
  scale_fill_manual(values = (brewer.pal(12, "Paired")[c(5,7,6,8)]))+
  scale_color_manual(values = (brewer.pal(12, "Paired")[c(5,7,6,8)]))+
  theme_few()+ theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))+
  labs(x = NULL, y = "Nuclei aspect ratio")

ggarrange(p1, p2, p3, ncol = 3)
