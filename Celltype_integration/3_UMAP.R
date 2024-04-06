library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(openxlsx)

job_id <- "651"

sobj <- readRDS("../data/sobj_x9_Wang_eLife_mnn_651.rds")

# Assuming 'celltype' and 'best_color' are vectors with your cell types and colors
# First, we ensure both are unique to match one-to-one as much as possible
celltype <- unique(sobj$celltype)
best_color<- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941",
               "#006FA6","#A30059","#FFE4E1","#0000A6","#63FFAC",
               "#B79762","#004D43","#8FB0FF","#997D87","#5A0007",
               "#809693","#1B4400","#4FC601","#3B5DFF","#FF2F80",
               "#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9",
               "#B903AA","#DDEFFF","#7B4F4B","#A1C299","#0AA6D8",
               "#00A087FF","#4DBBD5FF","#E64B35FF","#3C5488FF","#F38400",
               "#A1CAF1", "#C2B280","#848482","#E68FAC", "#0067A5", 
               "#F99379", "#604E97","#F6A600", "#B3446C","#DCD300",
               "#882D17", "#8DB600","#654522", "#E25822", "#2B3D26",
               "#191970","#000080",
               "#6495ED","#1E90FF","#00BFFF","#00FFFF","#FF1493",
               "#FF00FF","#A020F0","#63B8FF","#008B8B","#54FF9F",
               "#00FF00","#76EE00","#FFF68F","Yellow1","Gold1",
               "DarkGoldenrod4","#FF6A6A","#FF8247","#FFA54F","#FF7F24",
               "#FF3030","#FFA500","#FF7F00","#FF7256","#FF6347",
               "#FF4500","#FF1493","#FF6EB4","#EE30A7","#8B008B")

set.seed(123) # Ensure reproducibility if random selection is needed
if(length(celltype) > length(best_color)) {
  # If there are more cell types than colors, randomly pick colors for extra cell types
  extra_needed <- length(celltype) - length(best_color)
  random_colors <- sample(best_color, extra_needed, replace=TRUE)
  # Combine the original colors with the randomly picked additional colors
  final_colors <- c(best_color, random_colors)
} else {
  # If we have enough colors, or more, just use the first 'n' colors where 'n' is the number of cell types
  final_colors <- best_color[1:length(celltype)]
}

names(final_colors) <- celltype

# 'final_colors' is now your named vector with cell types and their colors
celltype_colors <- final_colors

# Create the plots
plot1 <- UMAPPlot(sobj, group.by = "celltype", label = TRUE, raster=FALSE) +
  scale_colour_manual(values = celltype_colors) +
  coord_fixed() +
  theme(legend.position = "none")  # Hide the legend for plot1

plot2 <- UMAPPlot(sobj, group.by = "celltype", split.by = "dataset", label = TRUE, raster=FALSE) +
  scale_colour_manual(values = celltype_colors) +
  coord_fixed()

# Combine the plots and add title
reference_plot <- plot1 + plot2 + plot_layout(widths = c(1, 4))
reference_plot <- reference_plot + 
  plot_annotation(title = "celltype", 
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 24)))
ggsave(paste0("../plots/UMAP_x9_Wang_eLife_", job_id, ".pdf"), reference_plot, width = 35, height = 7)

md <- sobj@meta.data
write.xlsx(md, file = "../data/md_x9_Wang_eLife.xlsx")

