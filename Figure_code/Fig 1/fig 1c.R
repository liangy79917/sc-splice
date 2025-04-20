library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(clustree)
library(cowplot)
library(qs)
library(bbknnR)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
options(future.globals.maxSize = 1000 * 1024^3)
sce <- readRDS("sce.rds")

data1 <- as.numeric(table(sce@meta.data$sample_ID))

sorted_data1 <- sort(data1)

median_value <- median(sorted_data1)
ggplot(data = data.frame(Index = 1:length(sorted_data1), Count = sorted_data1), aes(x = Index, y = Count)) +
  geom_bar(stat = 'identity', fill = "steelblue") +  
  geom_hline(yintercept = median_value, color = "black", linetype = "dashed") +  
#  scale_y_continuous(limits = c(0, 2500)) + 
  theme_minimal(base_size = 15) +  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "Individual index", y = "Cell number")

ggsave("cell_count.pdf", width = 8, height = 5)

print(paste("Median:", median_value))
