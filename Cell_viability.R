library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(DescTools)
library(ggpubr)



Viab_df <- read.csv("Skole/UiB/MCD-model/RNA-seq/Cell_viability.csv")
Assay_control_df <- read.csv("Skole/UiB/MCD-model/RNA-seq/Cell_viability_assay_control.csv")



#LPS
LPS_dat <- subset(Viab_df, Viab_df$Experiment == "LPS")



LPS_dat$Condition <- factor(LPS_dat$Condition, 
                            levels = c("Ctrl", "LPS_20", "LPS_40", "LPS_60",
                                       "LPS_80", "LPS_100", "LPS_120"))



LPS_aov <- aov(RFU ~ Condition, LPS_dat)
summary(LPS_aov)
DunnettTest(RFU ~ Condition, LPS_dat)



p1 <- ggplot(LPS_dat, aes(x = factor(Concentration), y = RFU, fill = Concentration))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch =21)+ 
  stat_summary(fun=mean, geom="crossbar")+
  scale_fill_distiller(palette = "Greys", direction = 1)+
  labs(x = "Concentration (µg/mL)", y = "Fluorescence\n(RFU)", subtitle = "LPS")+
  theme_few()+ theme(legend.position = "none")+
  ylim(0, 45000)+
  annotate("text", label = "ANOVA: ns", x = 1.5, y = 43000)



#PP2
PP2_dat <- subset(Viab_df, Viab_df$Experiment == "PP2'")



PP2_dat$Condition <- factor(PP2_dat$Condition, 
                            levels = c("Ctrl", "PP2_15", "PP2_30", "PP2_45",
                                       "PP2_60", "PP2_75", "PP2_90"))



PP2_aov <- aov(RFU ~ Condition, PP2_dat)
summary(PP2_aov)
DunnettTest(RFU ~ Condition, PP2_dat)



p2 <- ggplot(PP2_dat, aes(x = factor(Concentration), y = RFU, fill = Concentration))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch =21)+ 
  stat_summary(fun=mean, geom="crossbar")+
  scale_fill_distiller(palette = "Greys", direction = 1)+
  labs(x = "Concentration (µM)", y = NULL, subtitle = "PP2")+
  theme_few()+ theme(legend.position = "none")+
  ylim(0, 45000)+
  annotate("text", label = "ANOVA: ***\n0 µM vs. 90 µM: *", x = 1.5, y = 43000)



#Assay control: cells treated with 50% methanol
Assay_control_df$Treatment <- factor(Assay_control_df$Treatment)



lm_AssayCtrl <- lm(Fluorescence ~ Treatment, Assay_control_df)
summary(lm_AssayCtrl)



p3 <- ggplot(Assay_control_df, aes(x = factor(Concentration), y = Fluorescence))+
  geom_dotplot(binaxis = "y", stackdir = "center", pch =21, fill = "gray20")+
  stat_summary(fun=mean, geom="crossbar")+
  stat_compare_means(comparisons = list(c("0", "50")), 
                     method = "t.test", label = "p.signif")+
  theme_few()+
  labs(x = "Methanol (%)", y = NULL, subtitle = "Assay control")+
  ylim(0, 45000)



#Combine panels
ggarrange(p1, p2, p3, ncol = 3, widths = c(1,1, 0.5), align = "hv")



#Summary tables: mean +- SD
summary_stats_LPS <- LPS_dat %>%
  group_by(Concentration) %>%
  summarize(
    mean_RFU = mean(RFU, na.rm = TRUE),
    sd_RFU = sd(RFU, na.rm = TRUE),
    .groups = "drop")



summary_stats_PP2 <- PP2_dat %>%
  group_by(Concentration) %>%
  summarize(
    mean_RFU = mean(RFU, na.rm = TRUE),
    sd_RFU = sd(RFU, na.rm = TRUE),
    .groups = "drop")