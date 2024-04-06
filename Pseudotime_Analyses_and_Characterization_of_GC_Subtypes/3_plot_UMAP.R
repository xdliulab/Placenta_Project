library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(viridis)

job_id <- "645"

sobj <- readRDS("../data/sobj_x9_GC_mnn.rds")

Idents(sobj) <- "celltype"
DimPlot(sobj, label = TRUE)

celltype <- unique(sobj$celltype)

color <- read.csv("../../2024-03-17_627-CellChat-x9-snRNA-2GCsubtype-integrate/data/color_celltype_240318.csv", row.names = 1)
setdiff(color$order, celltype)
setdiff(celltype, color$order)

color <- color[color$order %in% celltype, ]

# Set the levels of the identity classes in your Seurat object
levels(sobj) <- color$order
color_vector <- setNames(color$Color, color$order)
p <- DimPlot(sobj, label = TRUE, cols = color_vector) + coord_fixed()
ggsave(filename = paste0("../plots/Trajecotry_GC_", job_id, ".pdf"), plot = p, 
       device = "pdf", width = 12, height = 10)

p <- DimPlot(sobj, label = FALSE, cols = color_vector) + 
  coord_fixed() +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(), # Explicitly removing axis lines
    legend.position = "none",
    panel.background = element_rect(fill = "transparent",color = NA), # Set panel background to transparent
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill='transparent', color=NA), # Set plot background to transparent
    legend.background = element_rect(fill='transparent'), # Ensure legend background is transparent
    legend.box.background = element_rect(fill='transparent') # Ensure legend box background is transparent
  )



ggsave(filename = paste0("../plots/Trajecotry_GC_nude_", job_id, ".png"), # Change the extension to .png
       plot = p, 
       device = "png", # Specify the device as PNG
       dpi = 300, # Set the resolution to 300 PPI
       width = 12, height = 10,
       bg = "transparent") # Ensure the background is transparent


# stage
Idents(sobj) <- "stage"
DimPlot(sobj, label = TRUE)
stage <- unique(sobj$stage)

colors <- c("#FDFFE1", "#E5F6D5", "#C3E7C7", "#94D2BA", "#58BAB1", "#469CAC", "#3377A1", "#244A88", "#24195B")
color_vector <- setNames(colors, levels(stage))
p <- DimPlot(sobj, label = TRUE, cols = color_vector) + coord_fixed()
ggsave(filename = paste0("../plots/Trajectory_GC_timepoints_", job_id, ".pdf"), plot = p, device = "pdf", 
       width = 12, height = 10)

# split stage
celltype <- unique(sobj$celltype)

color <- read.csv("../../2024-03-17_627-CellChat-x9-snRNA-2GCsubtype-integrate/data/color_celltype_240318.csv", row.names = 1)
setdiff(color$order, celltype)
setdiff(celltype, color$order)

color <- color[color$order %in% celltype, ]

# Set the levels of the identity classes in your Seurat object

Idents(sobj) <- "celltype"
levels(sobj) <- color$order
color_vector <- setNames(color$Color, color$order)

plot <- DimPlot(sobj, reduction = "umap", group.by = "celltype", split.by = "stage", cols = color_vector) + coord_fixed()

# Save the plot
ggsave(filename = paste0("../plots/Trajectory_GC_split_stage_", job_id, ".pdf"), 
       plot = plot, device = "pdf", width = 24, height = 4)
