

rm(list = ls())
gc()

library(tidyverse)
library(data.table)
library(fs)
library(stringr)
library(dplyr)
library(ggpmisc)
plot_data_1 <- read.table("plot_data_1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")

plot_data_2 <- plot_data_1 %>%
  mutate(sqtl = ifelse(sqtl >= 3, ">=3", as.character(sqtl))) %>%
  mutate(sqtl = factor(sqtl, levels = c("1", "2", ">=3")))

plot_data_3 <- plot_data_2 %>%
  group_by(cell_type) %>%
  summarize(sgene = n())

cell_type_order <- plot_data_3 %>%
  arrange(desc(sgene)) %>%
  select(cell_type) %>%
  unlist() %>%
  unname()
cell_type_order <- plot_data_3 %>%
  arrange(desc(sgene)) %>%
  select(cell_type) %>%
  unlist() %>%
  unname()
custom_sqtl_color <- c(
  "1" = "#FDE725", "2" = "#7AD151",
  ">3" = "#22A884"
)    
norm_const <- max(plot_data_3$sgene) * 1.10
out_dir <- "pdf/"
dir.create(out_dir, showWarnings = FALSE)
setwd("pdf/")

# Set custom order for cell types
cell_type_order <- c("BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg", 
                     "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT", 
                     "NKDim", "NKBright", "MonocM", "NoClaM")

# Create the plot
p1 <- ggplot() +
  # Stacked bar plot with custom color scheme
  geom_bar(
    data = plot_data_2,
    aes(x = cell_type, fill = sqtl),
    position = "fill",
    color = "black",
    linewidth = 0.3,
    width = 0.85
  ) +
  # Add points with labels
  geom_point(
    data = plot_data_3,
    aes(x = cell_type, y = sgene / norm_const),
    color = "#6A3D9A",
    size = 7,
    shape = 18
  ) +
  geom_text(
    data = plot_data_3,
    aes(x = cell_type, y = sgene / norm_const, label = sgene),
    color = "#6A3D9A",
    vjust = -1.8,
    size = 3.2,
    fontface = "bold"
  ) +
  # Y-axis settings
  scale_y_continuous(
    name = "Proportion of sGenes",
    expand = c(0, 0),
    labels = scales::percent_format(),
    sec.axis = sec_axis(
      trans = ~ . * norm_const,
      name = "Number of sGenes",
      breaks = scales::pretty_breaks(n = 6)
    )
  ) +
  scale_x_discrete(
    limits = cell_type_order  # Apply custom order to x-axis
  ) +
  # Color scheme and theme settings
  scale_fill_manual(
    values = c("#FDE725", "#7AD151", "#22A884"),  
    name = "sQTL/sGene",
    labels = c("1", "2", "≥3")
  ) +
  labs(
    title = "sQTL Distribution Across Cell Types",
    subtitle = "Stacked bars show proportion distribution, diamond points indicate total counts",
    caption = "Data Source: sQTL Analysis Pipeline"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "black", margin = margin(b = 15)),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      color = "black",
      size = 10
    ),
    axis.text.y = element_text(color = "black", size = 20),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "cm"),
    axis.line = element_line(color = "black"),
    plot.caption = element_text(color = "black", hjust = 1, margin = margin(t = 10)),
    axis.line.y.right = element_line(color = "#6A3D9A"),
    axis.ticks.y.right = element_line(color = "#6A3D9A"),
    axis.text.y.right = element_text(color = "#6A3D9A", size = 20),
    axis.title.y.right = element_text(color = "#6A3D9A", margin = margin(l = 10))
  ) +

  annotate(
    "segment",
    x = -Inf, xend = Inf,
    y = -Inf, yend = -Inf,
    color = "gray60", size = 0.8
  )

ggsave("sQTL.pdf", p1, width = 16, height = 10)

