library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)

job_id <- "645"

sobj <- readRDS("../data/sobj_x9_GC_mnn.rds")
DimPlot(sobj, group.by = "celltype")

cds <- as.cell_data_set(sobj)
cds <- cluster_cells(cds, resolution=1e-3)
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 10]))

plot <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "lightgrey",
           trajectory_graph_segment_size = 0) + coord_fixed()
ggsave(paste0("../plots/Trajectory_GC_pseudotime_", job_id, ".pdf"), plot, width = 12, height = 10, units = "in")



