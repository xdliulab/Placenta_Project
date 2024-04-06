library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(loupeR)

job_id <- "651"

# # Command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# 
# # Integration method and identifier
# integration_method <- args[1]
# identifier <- args[2]

integration_method <- "FastMNNIntegration"
identifier <- "mnn"
reduction_name <- paste('integrated', identifier, sep = '.')
sobj <- readRDS("../data/sobj_x9_Wang_eLife_before_integrate_mnn.rds")

# File paths
outRDS <- paste0("../data/sobj_x9_Wang_eLife_", identifier, "_", job_id, ".rds")
outPDF_QC <- paste0("../data/sobj_x9_Wang_eLife_", identifier, "_QC.pdf")
outPDF <- paste0("../data/sobj_x9_Wang_eLife_", identifier, ".pdf")
outLoupe <- paste0("../output/loupe_x9_Wang_eLife_", identifier, "_", job_id)

# Print the parameters
cat("Integration Method:", integration_method, "\n")
cat("Identifier:", identifier, "\n")
cat("Reduction Name:", reduction_name, "\n")

# Print file paths
cat("Output RDS File:", outRDS, "\n")
cat("Output QC PDF File:", outPDF_QC, "\n")
cat("Output PDF File:", outPDF, "\n")
cat("Output Loupe File:", outLoupe, "\n")

# Perform integration dynamically
sobj <- IntegrateLayers(
  object = sobj, method = integration_method,
  new.reduction = reduction_name,
  verbose = TRUE)

# Re-join layers after integration
sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])

# Further analysis steps
sobj <- FindNeighbors(sobj, reduction = reduction_name, dims = 1:30)
sobj <- FindClusters(sobj, resolution = seq(0.1, 0.1, 0.1))
sobj <- RunUMAP(sobj, dims = 1:30, reduction = reduction_name)

# Save plots
plot <- DimPlot(sobj, reduction = "umap", split.by = "dataset")
ggsave(outPDF_QC, plot, width = 7, height = 7, units = "in")

plot <- DimPlot(sobj, reduction = "umap", group.by = "celltype")
ggsave(outPDF, plot, width = 7, height = 7, units = "in")

# clusters <- select_clusters(sobj)
# projections <- select_projections(sobj)
# count_mat <- LayerData(sobj, assay = "RNA", layer = "counts")
# create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = outLoupe, force = TRUE)

# Save the integrated Seurat object
saveRDS(sobj, file = outRDS)

# Idents(sobj) <- "celltype"
# DimPlot(sobj, label = TRUE)
# celltype <- unique(sobj$celltype)
# 
# color <- read.csv("../../2024-01-24_495-x9-allcelltype-reclustering/data/color_celltype_240125.csv", row.names = 1)
# setdiff(color$order, celltype)
# setdiff(celltype, color$order)
# 
# # Set the levels of the identity classes in your Seurat object
# levels(sobj) <- color$order
# color_vector <- setNames(color$Color, color$order)
# p <- DimPlot(sobj, label = TRUE, cols = color_vector)
# ggsave(filename = "../plots/UMAP_allcelltype_v4.pdf", plot = p, 
#        device = "pdf", width = 14, height = 10)
