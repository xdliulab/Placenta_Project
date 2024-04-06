library(Seurat)
library(tidyverse)
library(patchwork)

sobj_x9 <- readRDS("../../2024-03-15_620-x9-snRNA-allcelltype-annotation-v4/data/sobj_x9_allcelltype_mnn_v4.rds")
sobj_Wang <- readRDS("../../2023-08-23_cell_subtype/data/2023-08-25_sobj_Wang.rds")
sobj_eLife <- readRDS("../../2023-08-23_cell_subtype/data/2023-08-25_sobj_eLife.rds")

DefaultAssay(sobj_x9) <- "RNA"
DefaultAssay(sobj_Wang) <- "RNA"
DefaultAssay(sobj_eLife) <- "RNA"

md_x9 <- sobj_x9@meta.data
md_Wang <- sobj_Wang@meta.data
md_eLife <- sobj_eLife@meta.data

desired_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "stage", "celltype")
sobj_x9@meta.data <- sobj_x9@meta.data[, desired_columns]
sobj_Wang@meta.data <- sobj_Wang@meta.data[, desired_columns]
sobj_eLife@meta.data <- sobj_eLife@meta.data[, desired_columns]

sobj_x9$orig.ident <- "x9"
sobj_Wang$orig.ident <- "Wang"
sobj_eLife$orig.ident <- "eLife"

sobj_x9$sample <- c(
  rep("x9_E9.5", sum(grepl("^E9.5", sobj_x9$stage))),
  rep("x9_E10.5", sum(grepl("^E10.5", sobj_x9$stage))),
  rep("x9_E11.5", sum(grepl("^E11.5", sobj_x9$stage))),
  rep("x9_E12.5", sum(grepl("^E12.5", sobj_x9$stage))),
  rep("x9_E13.5", sum(grepl("^E13.5", sobj_x9$stage))),
  rep("x9_E14.5", sum(grepl("^E14.5", sobj_x9$stage))),
  rep("x9_E15.5", sum(grepl("^E15.5", sobj_x9$stage))),
  rep("x9_E16.5", sum(grepl("^E16.5", sobj_x9$stage))),
  rep("x9_E18.5", sum(grepl("^E18.5", sobj_x9$stage)))
)

sobj_Wang$sample <- c(
  rep("Wang_E9.5", sum(grepl("-3$", colnames(sobj_Wang)))),
  rep("Wang_E10.5", sum(grepl("-5$", colnames(sobj_Wang)))),
  rep("Wang_E10.5h", sum(grepl("-6$", colnames(sobj_Wang)))),
  rep("Wang_E13.5", sum(grepl("-10$", colnames(sobj_Wang)))),
  rep("Wang_E14.5", sum(grepl("-11$", colnames(sobj_Wang))))
)

sobj_eLife$sample <- c(
  rep("eLife_E9.5", sum(grepl("^E9.5", sobj_eLife$stage))),
  rep("eLife_E10.5", sum(grepl("^E10.5", sobj_eLife$stage))),
  rep("eLife_E12.5", sum(grepl("^E12.5", sobj_eLife$stage))),
  rep("eLife_E14.5", sum(grepl("^E14.5", sobj_eLife$stage)))
)

sobj_x9[["RNA3"]] <- as(object = sobj_x9[["RNA"]], Class = "Assay")
DefaultAssay(sobj_x9) <- "RNA3"
sobj_x9[["RNA"]] <- NULL
sobj_x9 <- RenameAssays(object = sobj_x9, RNA3 = 'RNA')
sobj_x9$nCount_RNA3 <- NULL
sobj_x9$nFeature_RNA3 <- NULL

sobj_x9$dataset <- "x9"
sobj_Wang$dataset <- "Wang"
sobj_eLife$dataset <- "eLife"

sobj <- merge(sobj_x9, list(sobj_Wang, sobj_eLife))
md <- sobj@meta.data

sobj <- CreateSeuratObject(counts = sobj@assays[["RNA"]]@counts, meta.data = sobj@meta.data)

sobj
sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0)

sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$sample)
sobj

# run standard anlaysis workflow
sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj)

md <- sobj@meta.data
saveRDS(sobj, file = "../data/sobj_x9_Wang_eLife_before_integrate_mnn.rds")

# # Split each Seurat object by the "sample" column in metadata
# list_BGI <- SplitObject(sobj_x9, split.by = "sample")
# list_Deng <- SplitObject(sobj_Deng, split.by = "sample")
# list_eLife <- SplitObject(sobj_eLife, split.by = "sample")
# list_Wang <- SplitObject(sobj_Wang, split.by = "sample")
# 
# ifnb.list <- c(list_BGI, list_Deng, list_eLife, list_Wang)
# 
# rm(sobj_x9, sobj_Deng, sobj_eLife, sobj_Wang)
# gc() # Clear unused memory space
# rm(list_BGI, list_Deng, list_eLife, list_Wang)
# gc() # Clear unused memory space
# 
# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# # select features that are repeatedly variable across datasets for integration run PCA on each
# # dataset using these features
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#   x <- ScaleData(x, features = features, verbose = FALSE)
#   x <- RunPCA(x, features = features, verbose = FALSE)
# })
# 
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
# 
# rm(ifnb.list)
# gc() # Clear unused memory space
# 
# # this command creates an 'integrated' data assay
# immune.combined <- IntegrateData(anchorset = immune.anchors)
# 
# # Remove objects to save memory
# rm(immune.anchors)
# gc() # Clear unused memory space
# 
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
# DefaultAssay(immune.combined) <- "integrated"
# 
# # Run the standard workflow for visualization and clustering
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
# immune.combined <- FindClusters(immune.combined, resolution = seq(0.5, 3, 0.1))
# 
# saveRDS(immune.combined, file = "../data/2023-08-25_sobj_All4.rds")
# write.csv(immune.combined@meta.data, file = "../data/2023-08-25_md_All4.csv")
# 
# immune.combined <- readRDS("../data/2023-08-25_sobj_All4.rds")
# 
# # List to store individual plots
# plot_list <- list()
# # Iterate through the resolutions
# for(res in seq(0.5, 1, by=0.1)) {
#   # Define the name of the metadata column containing the clustering for this resolution
#   resolution_column <- paste0("integrated_snn_res.", res)
#   # Plot the UMAP for the current resolution
#   p <- DimPlot(immune.combined, group.by = resolution_column, raster=FALSE) + ggtitle(paste("Resolution", res))
#   # Add the plot to the list
#   plot_list[[length(plot_list) + 1]] <- p
# }
# # Combine the plots using patchwork
# combined_plot <- wrap_plots(plot_list)
# # Show the combined plot
# combined_plot
# # Save the combined plot as a PDF with the specified width and height
# ggsave("../plots/combined_plot_All4.pdf", combined_plot, width = 21, height = 14)
# 
# DimPlot(immune.combined, group.by = "orig.ident", raster = FALSE)
# DimPlot(immune.combined, group.by = "orig.ident", split.by = "orig.ident", raster = FALSE)
# 
# 
# 
# sobj <- readRDS("../data/2023-08-25_sobj_All4.rds")
# md <- sobj@meta.data
# DimPlot(sobj, group.by = "integrated_snn_res.3", raster = FALSE, label = TRUE)

