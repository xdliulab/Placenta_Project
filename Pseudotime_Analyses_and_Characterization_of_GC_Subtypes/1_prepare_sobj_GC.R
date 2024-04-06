library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)

sobj <- readRDS(paste0("../../2024-03-17_627-CellChat-x9-snRNA-2GCsubtype-integrate/output/sobj_snRNA_x9_allcelltype_2GCsubtypes_627.rds"))
md <- sobj@meta.data
selected_cluster <- c("JZP", "GC precursor", "GC-Prl7b1", "GC-Aldh1a3")
setdiff(selected_cluster, sobj$celltype)
sobj <- subset(sobj, subset = celltype %in% selected_cluster)
# sobj$celltype <- sobj$celltype
Idents(sobj) <- "celltype"
DimPlot(sobj, label = TRUE)

counts <- LayerData(sobj, assay = "RNA", layer = "counts")
sobj <- CreateSeuratObject(counts = counts, meta.data = sobj@meta.data)
sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$stage)
sobj

# run standard anlaysis workflow
sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj)

saveRDS(sobj, file = "../data/sobj_GC_before_integrate.rds")
