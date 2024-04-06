library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(loupeR)
library(viridis)

# # Command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# 
# # Integration method and identifier
# integration_method <- args[1]
# identifier <- args[2]

integration_method <- "FastMNNIntegration"
identifier <- "mnn"
reduction_name <- paste('integrated', identifier, sep = '.')
sobj <- readRDS("../data/sobj_GC_before_integrate.rds")

# File paths
outRDS <- paste0("../data/sobj_x9_GC_", identifier, ".rds")
outPDF_QC <- paste0("../data/sobj_x9_GC_", identifier, "_QC.pdf")
outPDF_stage <- paste0("../data/sobj_x9_GC_", identifier, "_stage.pdf")
outPDF <- paste0("../data/sobj_x9_GC_", identifier, ".pdf")
outLoupe <- paste0("../output/louper_x9_GC_", identifier, "")

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
plot <- DimPlot(sobj, reduction = "umap", split.by = "stage")
ggsave(outPDF_QC, plot, width = 7, height = 7, units = "in")

plot <- DimPlot(sobj, reduction = "umap", group.by = "stage", cols = viridis::viridis(9))  # Adjust the number inside viridis() based on your number of stage
ggsave(outPDF_stage, plot, width = 7, height = 7, units = "in")

plot <- DimPlot(sobj, reduction = "umap", group.by = "celltype", label = TRUE)
ggsave(outPDF, plot, width = 7, height = 7, units = "in")

clusters <- select_clusters(sobj)
projections <- select_projections(sobj)
count_mat <- sobj@assays[["RNA"]]@layers[["counts"]]
colnames(count_mat) <- colnames(sobj)
row.names(count_mat) <- row.names(sobj)
create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = outLoupe, force = TRUE)

# Save the integrated Seurat object
saveRDS(sobj, file = outRDS)
