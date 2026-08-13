library(RColorBrewer)
library(emmeans)
library(dplyr)
library(RColorBrewer)
library(WebPower)
library(ggplot2)
library(ggthemes)


Perm_df <- read.csv(file = "C:/Users/SimenFemangerPetters/OneDrive - Trosvik Consulting/Dokumenter/Skole/UiB/MCD-model/RNA-seq/Permeability.csv")

Perm_df <- Perm_df %>%
  mutate(Condition = factor(Condition, levels = c("Ctrl", "LPS trt", "PP2 trt")),
         Time = factor(Time))

head(Perm_df, 12)
str(Perm_df)


#Two-way repeated measures ANOVA
anova_model <- aov(Papp ~ Condition * Time + Error(Chip/Time), data = Perm_df)
summary(anova_model)



#Power analysis
#1) Calculate the effect size
F_val <- summary(anova_model)$`Error: Chip:Time`[[1]][["F value"]][[2]]
DF_effect <- summary(anova_model)$`Error: Chip:Time`[[1]][["Df"]][[2]]
DF_residual <- summary(anova_model)$`Error: Chip:Time`[[1]][["Df"]][[3]]

eta_sq <- (F_val * DF_effect) / ((F_val * DF_effect) + DF_residual)
Cohens_F <- sqrt(eta_sq / (1 - eta_sq))

print(c(F_val, DF_effect, DF_residual, eta_sq, Cohens_F))

#Run the power analysis
Power <- WebPower::wp.rmanova(n = 6,
                              ng = 3,
                              nm = 3,
                              f = Cohens_F,
                              type = 2)
Power



#Pairwise comparisons with p-value adjustment
emm <- emmeans(anova_model, ~ Condition | Time)
pairs(emm, adjust = "bonferroni")



#Plot
dodge = 0.1

plot_df <- Perm_df %>%
  mutate(X = as.numeric(Time), 
         X = case_when(Condition == "Ctrl" ~ X - dodge, Condition == "LPS trt" ~ X, TRUE ~ X + dodge))

ggplot(plot_df, aes(x = X, y = Papp, fill = Condition, group = Chip))+
  geom_line(linewidth = 0.2, color = "gray55") +
  geom_point(pch = 21, size = 4)+
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(5,7,8)])+
  labs(x = "Treatment time (hrs)", y = expression(P[app] ~ (cm/s)))+
  scale_x_continuous(breaks = 1:3, labels = c("0\n(Day 6)", "24\n(Day 7)", "48\n(Day 8)"))+
  annotate("text", label = "Condition: **\nTime : ***\nCondition:Time: ***", x = 1.3, y = 2.3e-05, size = 3)+
  theme_few()
