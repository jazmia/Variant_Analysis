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

cbPalette <- c("#0072B2","#0D0887","#9a179b", "#5b02a3","#CC79A7")
 
#Load the data
RIF <- read_excel("Sequencing_freq_data.xlsx", sheet = "RIF")

RIF$MIC <- factor(RIF$MIC, levels = c("0X", "2X", "4X", "16X", "32X", "64X"))

RIF$Mutation <- factor(RIF$Mutation, 
                       levels = c("rpoB", "inhA", "ppsD", "ppsE", "mptC"))

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
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", face = "bold", size = 10),
    plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10, face = "italic"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

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
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", face = "bold", size = 10),
    plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size = 14),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10, face = "italic"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

print(p) 

ggsave(filename = paste0("./graphs/LZD.svg"), 
       plot = p, width = 85, height=62.5, units = 'mm')
