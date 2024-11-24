#RIF fabF Growth Dynamics

# NOTES - normalised to min value

#Load packages 
library(tidyverse)
library(reshape2)
library(ggplot2)
library(readxl)
library(growthcurver)
library(openxlsx)
library(stringr)
library(ggtext)
library(rstatix)
library(ggpubr)
library(svglite)

#### RIF MIC GC Analysis ####
#### Formating Data ####

# Loading Data

fabF_Data <- read_excel("20240913_RIF_fabF.xlsx")


# Reformatting tables

long_fabF_Data <- reshape2::melt(fabF_Data, id.vars = "Hour")

long_fabF_Data <- long_fabF_Data %>% rename(Well = variable)


RIF_Layout <- read_excel("Plate_Layout.xlsx", sheet = "RIF")


fabF_Compiled_Data <- merge(long_fabF_Data, RIF_Layout, by = c("Well"))


# Antibiotic Concentration as factor
fabF_Compiled_Data$`Antibiotic Concentration` <- 
  factor(fabF_Compiled_Data$`Antibiotic Concentration`,
         levels = rev(c(20,18,16,14,12,10,8,6,4,2, 0, "Sterility Control")))

fabF_Compiled_Data <- fabF_Compiled_Data %>%
  mutate(Sample = case_when(
    Sample == "E. coli K12" ~ "WT",
    Sample == "E. coli K12 ΔfabF" ~ "ΔfabF",
    TRUE ~ Sample  # Keep the original value if it doesn't match any of the above
  ))

#### RIF FabF ####
# Graph each well
p = ggplot(fabF_Compiled_Data, aes(x = Hour, y = value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`)

print(p)

ggsave("./graphs/240913_fabF_overview.jpg", width = 297, 
       height = 210, units = "mm")

# Remove 40 conc antibiotic blank A02,  C03, C05

fabF_filtered_data <- fabF_Compiled_Data %>% 
  filter(!(Well %in% c("A02", "C03", "C05"))) %>%
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


#### Graph each well ####
p = ggplot(fabF_data, aes(x = Hour, y = normalised_value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`)

print(p)

ggsave("./graphs/240913_fabF_no_outlier_overview.jpg", width = 297, 
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
  scale_fill_manual(values = c("WT" = "#00BFC4", "ΔfabF" = "#7CAE00")) +
  scale_color_manual(values = c("WT" = "#00BFC4", "ΔfabF" = "#7CAE00"))

print(plot)

ggsave(filename = paste0("./graphs/240913_RIF_fabF_MIC_growth_curves.jpg"), 
       plot = plot, width = 297, height=210, units = 'mm')


#### AUC Analysis ####

#### Loading the Data ####

fabF_Data <- read_excel("20240913_RIF_fabF.xlsx")

## Rename Hour column to Time

colnames(fabF_Data)[1] <- "Time" 

## Summarise Data with Growthcurver

#Summarize the growth curve data within each well of a microtitre plate ---
summ_fabF_Data <- SummarizeGrowthByPlate(fabF_Data, bg_correct = "min") %>%
  select(sample, auc_l, note)

# Loading Plate Layout

fabF_Layout <- read_excel("Plate_Layout.xlsx", sheet = "RIF") %>%
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
  filter(Sample == "E. coli K12" & `Antibiotic Concentration` == 16)

K12_16_imputation <- K12_16 %>%
  summarise(auc_l = mean(auc_l)) %>%
  mutate(
    sample = "C03",
    note = "",
    `Antibiotic Concentration` = "16",
    Sample = "E. coli K12",
    Replicate = "R2"
  )

K12_12 <- rm_missing_data %>%
  filter(Sample == "E. coli K12" & `Antibiotic Concentration` == 12)

K12_12_imputation <- K12_12 %>%
  summarise(auc_l = mean(auc_l)) %>%
  mutate(
    sample = "C05",
    note = "",
    `Antibiotic Concentration` = "12",
    Sample = "E. coli K12",
    Replicate = "R2"
  )

missing_data <- rbind(K12_16_imputation, K12_12_imputation)

All_AUC_fabF <- rbind(rm_missing_data, missing_data)

# Organising data

All_AUC_fabF$`Antibiotic Concentration` <- 
  factor(All_AUC_fabF$`Antibiotic Concentration`, 
                              levels = rev(c("20",
                                         "18",
                                         "16",
                                         "14",
                                         "12",
                                         "10",
                                         "8",
                                         "6",
                                         "4",
                                         "2",
                                         "0",
                                         "Sterility Control")))

All_AUC_fabF$Sample <- 
  factor(All_AUC_fabF$Sample, 
         levels = c("E. coli K12", "E. coli K12 ΔfabF"))

#### T-test ####
# Boxplots for normality

ggplot(All_AUC_fabF, aes(x = `Antibiotic Concentration`, y = auc_l, 
                         fill = Sample)) +
  geom_boxplot(position = position_dodge(width = 0.75)) + 
  geom_point(position = position_dodge(width = 0.75), size = 1.25) +
  labs(title = "Normality Check",
       x = "Strain",
       y = "Area Under the Curve (AUC)") +
  theme(title = ggtext::element_markdown(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 55, hjust = 1))

ggsave(filename = "./t_test_boxplots/RIF.jpeg", 
       , width = 297, height=210, units = 'mm')

# Only three replicates so decent normality considering

sink("./t_test_assumptions/RIF.txt")
# Checking normality with Shapiro_Wilk
print(All_AUC_fabF %>%
        group_by(Sample, `Antibiotic Concentration`) %>%  
        shapiro_test(auc_l), n = Inf)
sink()

print(ggqqplot(All_AUC_fabF, "auc_l", ggtheme = theme_bw()) +
        facet_grid(`Antibiotic Concentration` ~ Sample, labeller = "label_both"))
ggsave(filename = "./t_test_assumptions/qqplot_RIF.jpeg", 
       , width = 297, height=210, units = 'mm')

sink("./t_test_assumptions/variance_RIF.txt")

# Checking homogneity of variance with levene_test
print(All_AUC_fabF %>% 
        levene_test(auc_l ~ interaction(`Antibiotic Concentration`, Sample)),
      n = Inf)
sink()

# T-test
sample_data <- All_AUC_fabF

t_test <- sample_data %>%
  group_by(`Antibiotic Concentration`) %>%
  t_test(auc_l ~ Sample) %>%
  add_xy_position(x = "Sample") %>%
  mutate(p.signif = case_when(
    p > 0.05 ~ "ns",
    p <= 0.05 & p > 0.01 ~ "*",
    p <= 0.01 & p > 0.001 ~ "**",
    p <= 0.001 & p > 0.0001 ~ "***",
    p <= 0.0001 ~ "****"
  ))

sink("./t_test/RIF.txt")
print(t_test)
sink()

#### Stats graphs ####
sample_data <- sample_data %>%
  mutate(Sample = case_when(
    Sample == "E. coli K12" ~ "WT",
    Sample == "E. coli K12 ΔfabF" ~ "ΔfabF",
    TRUE ~ Sample  # Keep the original value if it doesn't match any of the above
  ))

p = ggplot(sample_data, aes(x = Sample, y = auc_l)) +
  geom_boxplot(position = position_dodge(width = 0.75), aes(fill = Sample)) + 
  geom_point(position = position_dodge(width = 0.75), size = 1.25) +
  facet_grid(~ `Antibiotic Concentration`)+
  labs(title = "*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*",
       x = "Sample", 
       y = "AUC",
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
  scale_fill_manual(values = c("WT" = "#00BFC4", "ΔfabF" = "#7CAE00")) +
  stat_pvalue_manual(t_test, 
                     label = "p.signif", 
                     tip.length = 0.03, size = 6, hide.ns = TRUE)
    
print(p)
    
ggsave(gsub("%","_", paste0("./stat_graphs/RIF.jpg")), width = 297, height = 210, units = "mm")

stat_sample_data <- sample_data %>%
  group_by(Sample, `Antibiotic Concentration`) %>%
  summarise(
    mean = mean(auc_l),
    sd = sd(auc_l),
    n = n(),
    se = sd / sqrt(n))
head(stat_sample_data)

stat_sample_data <- stat_sample_data %>% rename(Strain = Sample)
head(stat_sample_data)

plot <- ggplot(stat_sample_data, 
               aes(x = `Antibiotic Concentration`, 
                   y = mean, colour = Strain, fill = Strain, group = Strain)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, 
                      ymax = mean + se), 
                  alpha = 0.7,
                width = 0.2) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "Antibiotic Conc. (µg/mL)", 
       y = "Mean AOU",
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
  scale_fill_manual(values = c("WT" = "#00BFC4", "ΔfabF" = "#7CAE00")) +
  scale_color_manual(values = c("WT" = "#00BFC4", "ΔfabF" = "#7CAE00"))

print(plot)

ggsave(filename = paste0("./graphs/240913_RIF_fabF_MIC_growth_curves_AUC.svg"), 
       plot = plot, width = 297, height=210, units = 'mm')

#### Difference ####
mean_data <- sample_data %>%
  group_by(Sample, `Antibiotic Concentration`) %>%
  summarise(mean_auc_l = mean(auc_l),
            se = sd(auc_l)/n())

difference_data <- mean_data %>%
  group_by(`Antibiotic Concentration`) %>%
  summarise(
    WT_mean = mean_auc_l[Sample == "WT"],
    fabF_mean = mean_auc_l[Sample == "ΔfabF"],
    Difference = fabF_mean - WT_mean,
    WT_se = se[Sample == "WT"],
    fabF_se = se[Sample == "ΔfabF"],
    SE_Difference = sqrt(WT_se^2 + fabF_se^2)
  )


relevant_stats <- t_test %>%
  select(`Antibiotic Concentration`, p, p.signif)

difference_data <- left_join(difference_data, relevant_stats, 
                             by = c("Antibiotic Concentration")) %>%
  mutate(y.position = case_when(difference_data$Difference >0 ~
                                Difference + 10,
                                difference_data$Difference <0 ~
                                Difference - 15)) %>%
  mutate(p.signif = case_when(
    p.signif == "ns" ~ "",
    TRUE ~ p.signif
  ),
  y.position = y.position + 3.5)

p = ggplot(difference_data, aes(x = factor(`Antibiotic Concentration`), y = Difference)) +
  geom_col(fill = "steelblue") +  # Create column plot
  geom_errorbar(aes(ymin = Difference - SE_Difference, ymax = Difference + SE_Difference), 
                width = 0.2, size = 0.7, color = "black") +
  geom_text(aes(label = p.signif, y = y.position), vjust = -0.5) +  # Add significance annotations
  geom_segment(aes(x = 7.7, xend = 8.3, y = 50, yend = 50), size = 0.5) +
  geom_segment(aes(x = 6.7, xend = 7.3, y = 97, yend = 97), size = 0.5) +
  geom_segment(aes(x = 5.7, xend = 6.3, y = 100, yend = 100), size = 0.5) + 
  geom_segment(aes(x = 4.7, xend = 5.3, y = 145, yend = 145), size = 0.5) + 
  geom_segment(aes(x = 3.7, xend = 4.3, y = 85, yend = 85), size = 0.5) + 
  labs(
    title = "Difference Between ΔfabF and WT Across Antibiotic Concentrations",
    x = "Antibiotic Concentration",
    y = "AUC Difference (ΔfabF - WT)"
  ) +
  theme(axis.text.x = element_text(size = 20),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.y = element_text(size = 20),  # Increase y-axis text size
        legend.title = element_text(size = 18),  # Larger legend title
        legend.text = element_text(size = 20),
        title = ggtext::element_markdown())

print(p)

ggsave("./difference_graphs/RIF.jpg", width = 297, height = 210, units = "mm")

#### Growth and Respiration Curves ####
#### Formating Data ####

# Loading Data

data <- read_excel("20241017_fabF_biolog.xlsx", sheet = "Data")

layout <- read_excel("20241017_fabF_biolog.xlsx", sheet = "Layout")

# Reformatting tables

long_data <- reshape2::melt(data, id.vars = "Hour")
long_data <- long_data %>% rename(Well = variable)

# Merging data and layout

Compiled_data <- merge(long_data, layout, by = c("Well"))

#### Overview ####

Compiled_data <- Compiled_data %>% 
  filter(`Antibiotic Concentration` != "Empty")

Compiled_data$`Antibiotic Concentration` <- factor(Compiled_data$`Antibiotic Concentration`,
                                                   levels = c("Sterility Control", "0", "5", "10","20"))

# View graph of each well

sample_data <- Compiled_data 

p = ggplot(sample_data, aes(x = Hour, y = value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`) +
  labs(y = "Strain", x = "Time")

print(p)

ggsave(gsub("%","_", paste0("./overview/sample_overview.jpg")), width = 297, height = 210, units = "mm")

sample_data <- Compiled_data %>%
  filter(`Fitness Measure`== "Growth")

p = ggplot(sample_data, aes(x = Hour, y = value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`) +
  labs(y = "Strain", x = "Time")

print(p)

ggsave(gsub("%","_", paste0("./overview/growth_overview.jpg")), width = 297, height = 210, units = "mm")
sample_data <- Compiled_data %>%
  filter(`Fitness Measure`== "Respiration")

p = ggplot(sample_data, aes(x = Hour, y = value, color = Well)) +
  geom_line() +
  facet_grid(Sample ~ `Antibiotic Concentration`) +
  labs(y = "Strain", x = "Time")

print(p)

ggsave(gsub("%","_", paste0("./overview/respiration_overview.jpg")), width = 297, height = 210, units = "mm")

## Filter out ##

filtered_data <- Compiled_data %>%
  filter(!(Well %in% c("A01", "C03", "D04", "C05")))

#### Normalisation ####
min_values <- filtered_data %>%
  group_by(Well, Sample, `Antibiotic Concentration`) %>%
  summarise(min_value = min(value)) %>%
  ungroup()

# Merge the minimum values back into the original data
normalised_data <- merge(filtered_data, min_values, 
                         by = c("Well", "Sample", "Antibiotic Concentration"))

# Normalize the value column by subtracting the minimum value
normalised_data$normalised_value <- normalised_data$value - normalised_data$min_value

# Remove the min_value column if it's no longer needed
normalised_data$min_value <- NULL

# Update fabF_data with the normalized data
data <- normalised_data %>%
  filter(Sample != "Blank")

#### Calculate  statistics ####

stat_data <- data %>%
  filter(`Antibiotic Concentration`!= "Sterility Control") %>%
  group_by(Hour, Sample, `Antibiotic Concentration`, `Fitness Measure`) %>%
  summarise(
    mean = mean(normalised_value),
    sd = sd(normalised_value),
    n = n(),
    se = sd / sqrt(n))
head(stat_data)

stat_data <- stat_data %>% rename(Strain = Sample)
head(stat_data)

#### Graphing Growth Curves ####

stat_data <- stat_data %>%
  mutate(Strain = case_when(
    Strain == "E. coli K12" ~ "WT",
    Strain == "E. coli K12 ΔfabF" ~ "ΔfabF",
    TRUE ~ Strain  # Keep the original value if it doesn't match any of the above
  ))

plot <- ggplot(stat_data, 
               aes(x = Hour, 
                   y = mean, 
                   color = `Antibiotic Concentration`,
                   fill = `Antibiotic Concentration`)) +
  geom_line() +
  facet_grid(`Fitness Measure`~ Strain, scales = "free_y") +
  geom_ribbon(aes(ymin = mean - se, 
                  ymax = mean + se), 
              alpha = 0.4, colour = NA) +
  labs(title = paste0("*E. coli* K12 WT vs *E. coli* K12 Δ*fabF*"),
       x = "Time (h)", 
       y = "Mean AOU",
       color = "Antibiotic Concentration") +
  theme(axis.text.x = element_text(size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.y = element_text(size = 20),  # Increase y-axis text size
        legend.title = element_text(size = 18),  # Larger legend title
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        title = ggtext::element_markdown())

print(plot)

ggsave(filename = paste0("./fabF_growth_respiration.jpg"), 
       plot = plot, width = 297, height=210, units = 'mm')
