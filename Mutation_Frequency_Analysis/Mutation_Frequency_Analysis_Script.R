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
library(viridis)

cbPalette <- c("#0D0887", "#9a179b", "#cb4679", "#ed7953", "#5b02a3")

#Load the data
RIF <- read_excel("Sequencing_freq_data.xlsx", sheet = "RIF")

RIF$MIC <- factor(RIF$MIC, levels = c("0X", "2X", "4X", "16X", "32X", "64X"))

RIF$Mutation <- factor(RIF$Mutation, 
                       levels = c("rpoB", "inhA", "ppsD", "ppsE", "Rv2181"))

LZD <- read_excel("Sequencing_freq_data.xlsx", sheet = "LZD")

LZD$MIC <- factor(LZD$MIC, levels = c("0X", "0.25X", "0.5X", "1X", "2X", "4X"))

LZD$Mutation <- factor(LZD$Mutation, 
                       levels = c("rplC", "glpK", "echA12"))

#Graph

p <- ggplot(RIF, aes(x = MIC, y = Frequency, 
                     colour = Mutation, shape = Mutation, group = Mutation)) +
  geom_line() +
  geom_point(size = 1.5) +
  labs(title = paste0("Rifampicin Directed Evolution"),
       x = "Relative MIC", 
       y = "Allele Frequency",
       color = "Mutation") +
  scale_shape_manual(values = c(0, 1, 2, 5, 6, 3)) +
  scale_colour_manual(values = cbPalette) +
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 10),  # Increase x-axis title size
        axis.title.y = element_text(size = 10),  # Increase y-axis title size
        axis.text.y = element_text(size = 9),  # Increase y-axis text size
        legend.title = element_text(size = 10),  # Larger legend title
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        title = ggtext::element_markdown()) 

print(p)

ggsave(filename = paste0("./graphs/RIF.svg"), 
       plot = p, width = 85, height=62.5, units = 'mm')

p <- ggplot(LZD, aes(x = MIC, y = Frequency, 
                     colour = Mutation, shape = Mutation, group = Mutation)) +
  geom_line() +
  geom_point(size = 1.5) +
  labs(title = paste0("Linezolid Directed Evolution"),
       x = "Relative MIC", 
       y = "Allele Frequency",
       color = "Mutation") +
  scale_shape_manual(values = c(0, 1, 2)) +
  scale_colour_manual(values = cbPalette) +
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 10),  # Increase x-axis title size
        axis.title.y = element_text(size = 10),  # Increase y-axis title size
        axis.text.y = element_text(size = 9),  # Increase y-axis text size
        legend.title = element_text(size = 10),  # Larger legend title
        legend.text = element_text(size = 10, face = "italic"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        title = ggtext::element_markdown())

print(p) 

ggsave(filename = paste0("./graphs/LZD.svg"), 
       plot = p, width = 85, height=62.5, units = 'mm')
