#### ANALYSIS OF WT and FABF E. COLI RIF MIC

# Notes - 
# normalised to min
# percentage growth calculate individually with mean of growth control

#Load packages 
library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(growthcurver)
library(stringr)
library(ggtext)
library(ggpubr)
library(svglite)
library(rstatix)

#### RIF MIC GC Analysis ####
#### fabF MIC Replicate 1 ####
#### Formating Data ####

#### Formating Data ####

# Loading Data

fabF_Data <- read_excel("250326_fabF_RIF_data.xlsx", sheet = "240913_R1")


# Reformatting tables

long_fabF_Data <- reshape2::melt(fabF_Data, id.vars = "Hour")

long_fabF_Data <- long_fabF_Data %>% rename(Well = variable)


RIF_Layout <- read_excel("250326_Plate_Layout.xlsx", sheet = "fabF")


fabF_Compiled_Data <- merge(long_fabF_Data, RIF_Layout, by = c("Well"))


# Antibiotic Concentration as factor
fabF_Compiled_Data$`Antibiotic Concentration` <- 
  factor(fabF_Compiled_Data$`Antibiotic Concentration`,
         levels = rev(c(20,18,16,14,12,10,8,6,4,2, 0, "Sterility Control")))

#### RIF FabF ####
# Graph each well
p = ggplot(fabF_Compiled_Data, aes(x = Hour, y = value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`)

print(p)

ggsave("./graphs/240913_fabF_overview.jpg", width = 297, 
       height = 210, units = "mm")

# Remove antibiotic blank, C03, D02

fabF_filtered_data <- fabF_Compiled_Data %>% 
  filter(!(Well %in% c("C03", "C05"))) %>%
  filter(Sample != "Blank")

fabF_data <- fabF_filtered_data

#### Normalisation ####
min_values <- fabF_filtered_data %>%
  group_by(Well, Sample, `Antibiotic Concentration`) %>%
  summarise(min_value = min(value)) %>%
  ungroup()

# Merge the minimum values back into the original data
normalised_fabF_data <- merge(fabF_filtered_data, min_values, 
                              by = c("Well", "Sample", "Antibiotic Concentration"))

# Normalize the value column by subtracting the minimum value
normalised_fabF_data$normalised_value <- normalised_fabF_data$value - normalised_fabF_data$min_value

# Remove the min_value column if it's no longer needed
normalised_fabF_data$min_value <- NULL

# Update fabF_data with the normalized data
fabF_data <- normalised_fabF_data

R1_data <- fabF_data


#### Graph each well ####
p = ggplot(fabF_data, aes(x = Hour, y = normalised_value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`)

print(p)

ggsave("./graphs/240913_fabF_cleaned_overview.jpg", width = 297, 
       height = 210, units = "mm")


#### Calculate  statistics ####

stat_fabF_data <- fabF_data %>%
  filter(`Antibiotic Concentration`!= "Sterility Control") %>%
  group_by(Hour, Sample, `Antibiotic Concentration`) %>%
  summarise(
    mean = mean(normalised_value),
    sd = sd(normalised_value),
    n = n(),
    se = sd / sqrt(n))
head(stat_fabF_data)

stat_fabF_data <- stat_fabF_data %>% rename(Strain = Sample)
head(stat_fabF_data)

#### Graphing Growth Curves ####

plot <- ggplot(stat_fabF_data, 
               aes(x = Hour, 
                   y = mean, 
                   color = Strain,
                   fill = Strain)) +
  geom_line() +
  facet_grid( ~ `Antibiotic Concentration`) +
  geom_ribbon(aes(ymin = mean - se, 
                  ymax = mean + se), 
              alpha = 0.4, colour = NA) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "Time (h)", 
       y = "Mean AOU",
       color = "Strain") +
  theme(axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.y = element_text(size = 20),  # Increase y-axis text size
        legend.title = element_text(size = 18),  # Larger legend title
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  scale_color_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b"))

print(plot)

ggsave(filename = paste0("./graphs/240913_RIF_fabF_MIC_growth_curves.jpg"), 
       plot = plot, width = 297, height=210, units = 'mm')


#### AUC Analysis ####

#### Loading the Data ####

fabF_Data <- read_excel("250326_fabF_RIF_data.xlsx", sheet = "240913_R1")

## Rename Hour column to Time

colnames(fabF_Data)[1] <- "Time" 

## Summarise Data with Growthcurver

#Summarize the growth curve data within each well of a microtitre plate ---
summ_fabF_Data <- SummarizeGrowthByPlate(fabF_Data, bg_correct = "min") %>%
  select(sample, auc_l, note)

# Loading Plate Layout

fabF_Layout <- read_excel("250326_Plate_Layout.xlsx", sheet = "fabF") %>%
  rename(sample = Well)

## Merging layout and summary

Compiled_fabF <- merge(summ_fabF_Data, fabF_Layout, by = c("sample"))

All_AUC_fabF <- Compiled_fabF %>% 
  filter(Sample != "Blank") %>%
  filter(`Antibiotic Concentration` != "Sterility Control")

#### Dealing with missing data ####

rm_missing_data <- All_AUC_fabF %>%
  filter(!(sample %in% c("C03", "C05")))

K12_16 <- rm_missing_data %>%
  filter(Sample == "WT" & `Antibiotic Concentration` == 16)

K12_16_imputation <- K12_16 %>%
  summarise(auc_l = mean(auc_l)) %>%
  mutate(
    sample = "C03",
    note = "",
    `Antibiotic Concentration` = "16",
    Sample = "WT",
    Replicate = "R2"
  )

K12_12 <- rm_missing_data %>%
  filter(Sample == "WT" & `Antibiotic Concentration` == 12)

K12_12_imputation <- K12_12 %>%
  summarise(auc_l = mean(auc_l)) %>%
  mutate(
    sample = "C05",
    note = "",
    `Antibiotic Concentration` = "12",
    Sample = "WT",
    Replicate = "R2"
  )

missing_data <- rbind(K12_16_imputation, K12_12_imputation)

All_AUC_fabF <- rbind(rm_missing_data, missing_data)

# Organising data

All_AUC_fabF$Sample <- 
  factor(All_AUC_fabF$Sample, 
         levels = c("WT", "ΔfabF")) 

All_AUC_fabF <- All_AUC_fabF %>%
  rename(Strain = Sample) %>%
  filter(`Antibiotic Concentration` != "Sterility Control")

All_AUC_fabF$`Antibiotic Concentration` <- 
  factor(All_AUC_fabF$`Antibiotic Concentration`,
         levels = rev(c(20,18,16,14,12,10,8,6,4,2, 0, "Sterility Control")))

#### Analysis % growth relative to GC ####
strain_GC <- All_AUC_fabF %>%
  filter(`Antibiotic Concentration` == 0) %>%
  select(Strain, auc_l)

strain_GC_stats <- strain_GC %>%
  group_by(Strain) %>%
  summarise(growth_control_mean = mean(auc_l))

# Compute percentage growth AUC
percent_growth_fabF <- All_AUC_fabF %>%
  left_join(strain_GC_stats, by = "Strain") %>%
  mutate(percent_growth = (auc_l / growth_control_mean) * 100)

percent_growth_R1 <- percent_growth_fabF

#### T-test ####
# Boxplots for normality

ggplot(percent_growth_fabF, aes(x = `Antibiotic Concentration`, y = percent_growth, 
                                fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.75)) + 
  geom_point(position = position_dodge(width = 0.75), size = 1.25) +
  labs(title = "Normality Check",
       x = "Strain",
       y = "Area Under the Curve (AUC)") +
  theme(title = ggtext::element_markdown(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 55, hjust = 1))

ggsave(filename = "./stats/240913_fabF_percent_growth.jpeg", 
       width = 297, height=210, units = 'mm')

# Only three replicates so decent normality considering

sink("./stats/240913_fabF.txt")
# Checking normality with Shapiro_Wilk
print(percent_growth_fabF %>%
        group_by(Strain, `Antibiotic Concentration`) %>%  
        shapiro_test(percent_growth), n = Inf)
sink()

print(ggqqplot(percent_growth_fabF, "percent_growth", ggtheme = theme_bw()) +
        facet_grid(`Antibiotic Concentration` ~ Strain, labeller = "label_both"))
ggsave(filename = "./stats/240913_fabF_qqplot.jpeg", 
       , width = 297, height=210, units = 'mm')

sink("./stats/240913_fabF_variance.txt")

# Checking homogneity of variance with levene_test
print(percent_growth_fabF %>% 
        levene_test(percent_growth ~ interaction(`Antibiotic Concentration`, Strain)),
      n = Inf)
sink()

# T-test

t_test <- percent_growth_fabF %>%
  group_by(`Antibiotic Concentration`) %>%
  t_test(percent_growth ~ Strain) %>%
  add_xy_position(x = "Strain") %>%
  mutate(p.signif = case_when(
    p > 0.05 ~ "ns",
    p <= 0.05 & p > 0.01 ~ "*",
    p <= 0.01 & p > 0.001 ~ "**",
    p <= 0.001 & p > 0.0001 ~ "***",
    p <= 0.0001 ~ "****"
  ))

sink("./stats/240913_fabF_t_test.txt")
print(t_test)
sink()

# Stats graphs

p = ggplot(percent_growth_fabF, aes(x = Strain, y = percent_growth)) +
  geom_boxplot(position = position_dodge(width = 0.75), aes(fill = Strain)) + 
  geom_point(position = position_dodge(width = 0.75), size = 1.25) +
  facet_grid(~ `Antibiotic Concentration`)+
  labs(title = "*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*",
       x = "Sample", 
       y = "Percentage AUC Growth (normalised to no drug)",
       caption = "Statistics: t-test") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 38),
        strip.text = element_text(size = 26),  # Text size for facet labels
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 32),
        plot.caption = element_text(size = 20),
        axis.line = element_line(color = "black"),  
        panel.border = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  stat_pvalue_manual(t_test, 
                     label = "p.signif", 
                     tip.length = 0.03, size = 6, hide.ns = TRUE)

print(p)

ggsave(gsub("%","_", paste0("./graphs/240913_fabF_stats.jpg")), width = 297, height = 210, units = "mm")

# summary stats

summary_stats_data <- percent_growth_fabF %>%
  group_by(Strain, `Antibiotic Concentration`) %>%
  summarise(
    mean = mean(percent_growth),
    sd = sd(percent_growth),
    n = n(),
    se = sd / sqrt(n))
head(summary_stats_data)

plot <- ggplot(summary_stats_data, 
               aes(x = `Antibiotic Concentration`, 
                   y = mean, colour = Strain, fill = Strain, group = Strain)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                alpha = 0.7,
                width = 0.2) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "RIF (µg/mL)", 
       y = "Percentage AUC Growth (normalised to no drug)",
       color = "Strain") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 38),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 32),
        axis.line = element_line(color = "black"),  
        panel.border = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  scale_color_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b"))

print(plot)

ggsave(filename = paste0("./graphs/240913_fabF_mean.svg"), 
       plot = plot, width = 297, height=210, units = 'mm')
#### ---------------------------------------------------------------- #####
#### fabF MIC Replicate 2 ####

#### Formating Data ####

#### Formating Data ####

# Loading Data

fabF_Data <- read_excel("250326_fabF_RIF_data.xlsx", sheet = "250326_R2")


# Reformatting tables

long_fabF_Data <- reshape2::melt(fabF_Data, id.vars = "Hour")

long_fabF_Data <- long_fabF_Data %>% rename(Well = variable)


RIF_Layout <- read_excel("250326_Plate_Layout.xlsx", sheet = "fabF")


fabF_Compiled_Data <- merge(long_fabF_Data, RIF_Layout, by = c("Well"))


# Antibiotic Concentration as factor
fabF_Compiled_Data$`Antibiotic Concentration` <- 
  factor(fabF_Compiled_Data$`Antibiotic Concentration`,
         levels = rev(c(20,18,16,14,12,10,8,6,4,2, 0, "Sterility Control")))

#### RIF FabF ####
# Graph each well
p = ggplot(fabF_Compiled_Data, aes(x = Hour, y = value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`)

print(p)

ggsave("./graphs/250326_fabF_overview.jpg", width = 297, 
       height = 210, units = "mm")

# Remove antibiotic blank, C03, D02

fabF_filtered_data <- fabF_Compiled_Data %>% 
  filter(!(Well %in% c("C03", "D02"))) %>%
  filter(Sample != "Blank")

fabF_data <- fabF_filtered_data

#### Normalisation ####
min_values <- fabF_filtered_data %>%
  group_by(Well, Sample, `Antibiotic Concentration`) %>%
  summarise(min_value = min(value)) %>%
  ungroup()

# Merge the minimum values back into the original data
normalised_fabF_data <- merge(fabF_filtered_data, min_values, 
                              by = c("Well", "Sample", "Antibiotic Concentration"))

# Normalize the value column by subtracting the minimum value
normalised_fabF_data$normalised_value <- normalised_fabF_data$value - normalised_fabF_data$min_value

# Remove the min_value column if it's no longer needed
normalised_fabF_data$min_value <- NULL

# Update fabF_data with the normalized data
fabF_data <- normalised_fabF_data

R2_data <- fabF_data

#### Graph each well ####
p = ggplot(fabF_data, aes(x = Hour, y = normalised_value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`)

print(p)

ggsave("./graphs/250326_fabF_cleaned_overview.jpg", width = 297, 
       height = 210, units = "mm")


#### Calculate  statistics ####

stat_fabF_data <- fabF_data %>%
  filter(`Antibiotic Concentration`!= "Sterility Control") %>%
  group_by(Hour, Sample, `Antibiotic Concentration`) %>%
  summarise(
    mean = mean(normalised_value),
    sd = sd(normalised_value),
    n = n(),
    se = sd / sqrt(n))
head(stat_fabF_data)

stat_fabF_data <- stat_fabF_data %>% rename(Strain = Sample)
head(stat_fabF_data)

#### Graphing Growth Curves ####

plot <- ggplot(stat_fabF_data, 
               aes(x = Hour, 
                   y = mean, 
                   color = Strain,
                   fill = Strain)) +
  geom_line() +
  facet_grid( ~ `Antibiotic Concentration`) +
  geom_ribbon(aes(ymin = mean - se, 
                  ymax = mean + se), 
              alpha = 0.4, colour = NA) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "Time (h)", 
       y = "Mean AOU",
       color = "Strain") +
  theme(axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.y = element_text(size = 20),  # Increase y-axis text size
        legend.title = element_text(size = 18),  # Larger legend title
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  scale_color_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b"))

print(plot)

ggsave(filename = paste0("./graphs/250326_RIF_fabF_MIC_growth_curves.jpg"), 
       plot = plot, width = 297, height=210, units = 'mm')


#### AUC Analysis ####

#### Loading the Data ####

fabF_Data <- read_excel("250326_fabF_RIF_data.xlsx", sheet = "250326_R2")

## Rename Hour column to Time

colnames(fabF_Data)[1] <- "Time" 

## Summarise Data with Growthcurver

#Summarize the growth curve data within each well of a microtitre plate ---
summ_fabF_Data <- SummarizeGrowthByPlate(fabF_Data, bg_correct = "min") %>%
  select(sample, auc_l, note)

# Loading Plate Layout

fabF_Layout <- read_excel("250326_Plate_Layout.xlsx", sheet = "fabF") %>%
  rename(sample = Well)

## Merging layout and summary

Compiled_fabF <- merge(summ_fabF_Data, fabF_Layout, by = c("sample"))

All_AUC_fabF <- Compiled_fabF %>% 
  filter(Sample != "Blank") %>%
  filter(`Antibiotic Concentration` != "Sterility Control")

#### Dealing with missing data ####

rm_missing_data <- All_AUC_fabF %>%
  filter(!(sample %in% c("C03", "D02")))

K12_16 <- rm_missing_data %>%
  filter(Sample == "WT" & `Antibiotic Concentration` == 16)

K12_16_imputation <- K12_16 %>%
  summarise(auc_l = mean(auc_l)) %>%
  mutate(
    sample = "C03",
    note = "",
    `Antibiotic Concentration` = "16",
    Sample = "WT",
    Replicate = "R2"
  )

K12_18 <- rm_missing_data %>%
  filter(Sample == "WT" & `Antibiotic Concentration` == 18)

K12_18_imputation <- K12_18 %>%
  summarise(auc_l = mean(auc_l)) %>%
  mutate(
    sample = "D02",
    note = "",
    `Antibiotic Concentration` = "18",
    Sample = "WT",
    Replicate = "R3"
  )

missing_data <- rbind(K12_16_imputation, K12_18_imputation)

All_AUC_fabF <- rbind(rm_missing_data, missing_data)

# Organising data

All_AUC_fabF$Sample <- 
  factor(All_AUC_fabF$Sample, 
         levels = c("WT", "ΔfabF")) 

All_AUC_fabF <- All_AUC_fabF %>%
  rename(Strain = Sample) %>%
  filter(`Antibiotic Concentration` != "Sterility Control")

All_AUC_fabF$`Antibiotic Concentration` <- 
  factor(All_AUC_fabF$`Antibiotic Concentration`,
         levels = rev(c(20,18,16,14,12,10,8,6,4,2, 0, "Sterility Control")))

#### Analysis % growth relative to GC ####
strain_GC <- All_AUC_fabF %>%
  filter(`Antibiotic Concentration` == 0) %>%
  select(Strain, auc_l)

strain_GC_stats <- strain_GC %>%
  group_by(Strain) %>%
  summarise(growth_control_mean = mean(auc_l))

# Compute percentage growth AUC
percent_growth_fabF <- All_AUC_fabF %>%
  left_join(strain_GC_stats, by = "Strain") %>%
  mutate(percent_growth = (auc_l / growth_control_mean) * 100) 

percent_growth_R2 <- percent_growth_fabF

#### T-test ####
# Boxplots for normality

ggplot(percent_growth_fabF, aes(x = `Antibiotic Concentration`, y = percent_growth, 
                                fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.75)) + 
  geom_point(position = position_dodge(width = 0.75), size = 1.25) +
  labs(title = "Normality Check",
       x = "Strain",
       y = "Area Under the Curve (AUC)") +
  theme(title = ggtext::element_markdown(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 55, hjust = 1))

ggsave(filename = "./stats/250326_fabF_percent_growth.jpeg", 
       width = 297, height=210, units = 'mm')

# Only three replicates so decent normality considering

sink("./stats/250326_fabF.txt")
# Checking normality with Shapiro_Wilk
print(percent_growth_fabF %>%
        group_by(Strain, `Antibiotic Concentration`) %>%  
        shapiro_test(percent_growth), n = Inf)
sink()

print(ggqqplot(percent_growth_fabF, "percent_growth", ggtheme = theme_bw()) +
        facet_grid(`Antibiotic Concentration` ~ Strain, labeller = "label_both"))
ggsave(filename = "./stats/250326_fabF_qqplot.jpeg", 
       , width = 297, height=210, units = 'mm')

sink("./stats/250326_fabF_variance.txt")

# Checking homogneity of variance with levene_test
print(percent_growth_fabF %>% 
        levene_test(percent_growth ~ interaction(`Antibiotic Concentration`, Strain)),
      n = Inf)
sink()

# T-test

t_test <- percent_growth_fabF %>%
  group_by(`Antibiotic Concentration`) %>%
  t_test(percent_growth ~ Strain) %>%
  add_xy_position(x = "Strain") %>%
  mutate(p.signif = case_when(
    p > 0.05 ~ "ns",
    p <= 0.05 & p > 0.01 ~ "*",
    p <= 0.01 & p > 0.001 ~ "**",
    p <= 0.001 & p > 0.0001 ~ "***",
    p <= 0.0001 ~ "****"
  ))

sink("./stats/250326_fabF_t_test.txt")
print(t_test)
sink()

# Stats graphs

p = ggplot(percent_growth_fabF, aes(x = Strain, y = percent_growth)) +
  geom_boxplot(position = position_dodge(width = 0.75), aes(fill = Strain)) + 
  geom_point(position = position_dodge(width = 0.75), size = 1.25) +
  facet_grid(~ `Antibiotic Concentration`)+
  labs(title = "*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*",
       x = "Sample", 
       y = "Percentage AUC Growth (normalised to no drug)",
       caption = "Statistics: t-test") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 38),
        strip.text = element_text(size = 26),  # Text size for facet labels
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 32),
        plot.caption = element_text(size = 20),
        axis.line = element_line(color = "black"),  
        panel.border = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  stat_pvalue_manual(t_test, 
                     label = "p.signif", 
                     tip.length = 0.03, size = 6, hide.ns = TRUE)

print(p)

ggsave(gsub("%","_", paste0("./graphs/250326_fabF_stats.jpg")), width = 297, height = 210, units = "mm")

# summary stats

summary_stats_data <- percent_growth_fabF %>%
  group_by(Strain, `Antibiotic Concentration`) %>%
  summarise(
    mean = mean(percent_growth),
    sd = sd(percent_growth),
    n = n(),
    se = sd / sqrt(n))
head(summary_stats_data)

plot <- ggplot(summary_stats_data, 
               aes(x = `Antibiotic Concentration`, 
                   y = mean, colour = Strain, fill = Strain, group = Strain)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                alpha = 0.7,
                width = 0.2) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "RIF (µg/mL)", 
       y = "Percentage AUC Growth (normalised to no drug)",
       color = "Strain") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 38),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 32),
        axis.line = element_line(color = "black"),  
        panel.border = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  scale_color_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b"))

print(plot)

ggsave(filename = paste0("./graphs/250326_fabF_mean.svg"), 
       plot = plot, width = 297, height=210, units = 'mm')
#### ---------------------------------------------------------------- #####
#### Combined replicates ####
# Loading and Reformatting Data

R1_data <- R1_data %>%
  mutate(bio_rep = "R1")

R2_data <- R2_data %>%
  mutate(bio_rep = "R2")

combined_data <- rbind(R1_data, R2_data)

#### Graph replicates ####
p = ggplot(combined_data, aes(x = Hour, y = normalised_value, color = bio_rep)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`)

print(p)

ggsave("./graphs/combined_overview.jpg", width = 297, 
       height = 210, units = "mm")

#### Calculate  statistics ####

stat_combined_data <- combined_data %>%
  filter(`Antibiotic Concentration`!= "Sterility Control") %>%
  group_by(Hour, Sample, `Antibiotic Concentration`) %>%
  summarise(
    mean = mean(normalised_value),
    sd = sd(normalised_value),
    n = n(),
    se = sd / sqrt(n))
head(stat_combined_data)

stat_combined_data <- stat_combined_data %>% rename(Strain = Sample)
head(stat_combined_data)

#### Graphing Growth Curves ####

plot <- ggplot(stat_combined_data, 
               aes(x = Hour, 
                   y = mean, 
                   color = Strain,
                   fill = Strain)) +
  geom_line() +
  facet_grid( ~ `Antibiotic Concentration`) +
  geom_ribbon(aes(ymin = mean - se, 
                  ymax = mean + se), 
              alpha = 0.4, colour = NA) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "Time (h)", 
       y = "Mean AOU",
       color = "Strain") +
  theme(axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.y = element_text(size = 20),  # Increase y-axis text size
        legend.title = element_text(size = 18),  # Larger legend title
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  scale_color_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b"))

print(plot)

ggsave(filename = paste0("./graphs/combined_MIC_growth_curves.jpg"), 
       plot = plot, width = 297, height=210, units = 'mm')

#### AUC Analysis ####

#### Loading the Data ####

percent_growth_R1 <- percent_growth_R1 %>%
  mutate(bio_rep = "R1")

percent_growth_R2 <- percent_growth_R2 %>%
  mutate(bio_rep = "R2")

combined_percent_growth <- rbind(percent_growth_R1, percent_growth_R2)

#### Analysis % growth relative to GC ####

#### T-test ####
# Boxplots for normality

ggplot(combined_percent_growth, aes(x = `Antibiotic Concentration`, y = percent_growth, linetype = Strain)) +
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(aes(colour = Strain, shape = bio_rep), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), 
              size = 2, stroke = 1) +
  labs(title = "Normality Check",
       x = "Strain",
       y = "Area Under the Curve (AUC)") +
  theme(title = ggtext::element_markdown(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 55, hjust = 1))

ggsave(filename = "./stats/combined_percent_growth.jpeg", 
       width = 297, height=210, units = 'mm')

# Only three replicates so decent normality considering

sink("./stats/combined.txt")
# Checking normality with Shapiro_Wilk
print(combined_percent_growth %>%
        group_by(Strain, `Antibiotic Concentration`) %>%  
        shapiro_test(percent_growth), n = Inf)
sink()

print(ggqqplot(combined_percent_growth, "percent_growth", ggtheme = theme_bw()) +
        facet_grid(`Antibiotic Concentration` ~ Strain, labeller = "label_both"))
ggsave(filename = "./stats/combined_qqplot.jpeg", 
       , width = 297, height=210, units = 'mm')

sink("./stats/combined_variance.txt")

# Checking homogneity of variance with levene_test - not necessary for Welch's t-test
print(combined_percent_growth %>% 
        levene_test(percent_growth ~ interaction(`Antibiotic Concentration`, Strain)),
      n = Inf)
sink()

# T-test

t_test <- combined_percent_growth %>%
  group_by(`Antibiotic Concentration`) %>%
  t_test(percent_growth ~ Strain) %>%
  add_xy_position(x = "Strain") %>%
  mutate(p.signif = case_when(
    p > 0.05 ~ "ns",
    p <= 0.05 & p > 0.01 ~ "*",
    p <= 0.01 & p > 0.001 ~ "**",
    p <= 0.001 & p > 0.0001 ~ "***",
    p <= 0.0001 ~ "****"
  ))

sink("./stats/combined_t_test.txt")
print(t_test)
sink()

# Stats graphs

p = ggplot(combined_percent_growth, aes(x = Strain, y = percent_growth)) +
  geom_boxplot(position = position_dodge(width = 0.75), aes(fill = Strain)) + 
  geom_point(position = position_dodge(width = 0.75), size = 1.25) +
  facet_grid(~ `Antibiotic Concentration`)+
  labs(title = "*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*",
       x = "Sample", 
       y = "Percentage AUC Growth (normalised to no drug)",
       caption = "Statistics: t-test") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 38),
        strip.text = element_text(size = 26),  # Text size for facet labels
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 32),
        plot.caption = element_text(size = 20),
        axis.line = element_line(color = "black"),  
        panel.border = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  stat_pvalue_manual(t_test, 
                     label = "p.signif", 
                     tip.length = 0.03, size = 6, hide.ns = TRUE)

print(p)

ggsave(gsub("%","_", paste0("./graphs/combined_stats.svg")), width = 297, height = 210, units = "mm")

# summary stats

summary_stats_data <- combined_percent_growth %>%
  group_by(Strain, `Antibiotic Concentration`) %>%
  summarise(
    mean = mean(percent_growth),
    sd = sd(percent_growth),
    n = n(),
    se = sd / sqrt(n))
head(summary_stats_data)

plot <- ggplot(summary_stats_data, 
               aes(x = `Antibiotic Concentration`, 
                   y = mean, colour = Strain, fill = Strain, group = Strain,
                   shape = Strain)) +
  geom_line() +
  geom_point(size = 1.25) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se), 
                alpha = 0.7,
                width = 0.2) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "RIF (µg/mL)", 
       y = "Percentage AUC Growth (normalised to no drug)",
       color = "Strain") +
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 10),  # Increase x-axis title size
        axis.title.y = element_text(size = 10),  # Increase y-axis title size
        axis.text.y = element_text(size = 9),  # Increase y-axis text size
        legend.title = element_text(size = 10),  # Larger legend title
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        title = ggtext::element_markdown()) +
  scale_fill_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  scale_color_manual(values = c("WT" = "#3E3E3E", "ΔfabF" = "#9a179b")) +
  scale_shape_manual(values = c(1,2))+
  ylim(0, NA)

print(plot)

ggsave(filename = paste0("./graphs/combined_mean.svg"), 
       plot = plot, width = 85, height=85, units = 'mm')

