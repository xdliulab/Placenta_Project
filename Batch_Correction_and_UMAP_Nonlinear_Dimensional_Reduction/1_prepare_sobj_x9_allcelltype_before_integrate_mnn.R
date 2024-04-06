library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)

counts <- Read10X(data.dir = paste0("../../2024-01-22_450-cellranger-10x9-aggr/data/count/filtered_feature_bc_matrix/"))
barcodes_allcelltype <- read.csv("../../2024-01-24_493-x9-allcelltype-annotation/output/allcelltype_updated_v2.csv")
# all.equal(colnames(counts), barcodes_allcelltype$Barcode)
unique(barcodes_allcelltype$RNA_snn_res.0.1)
outRDS <- "../data/sobj_x9_allcelltype_before_integrate_mnn.rds"

filtered_df <- barcodes_allcelltype
# # Filter the dataframe
# filtered_df <- barcodes_allcelltype %>%
#   filter(RNA_snn_res.0.1 %in% c("Trophoblast-a",
#                                 "Trophoblast-b",
#                                 "Trophoblast-c",
#                                 "Trophoblast-d",
#                                 "Trophoblast-e",
#                                 "Trophoblast-f"))

filtered_counts <- counts[, filtered_df$Barcode]

# Initialize the Seurat object with the raw (non-normalized data).
sobj <- CreateSeuratObject(counts = filtered_counts)
sobj
sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0)

sobj$stages <- c(
  rep("E9.5", sum(grepl("-1$", colnames(sobj)))),
  rep("E10.5", sum(grepl("-2$", colnames(sobj)))),
  rep("E11.5", sum(grepl("-3$", colnames(sobj)))),
  rep("E12.5", sum(grepl("-4$", colnames(sobj)))),
  rep("E13.5", sum(grepl("-5$", colnames(sobj)))),
  rep("E14.5", sum(grepl("-6$", colnames(sobj)))),
  rep("E15.5", sum(grepl("-7$", colnames(sobj)))),
  rep("E16.5", sum(grepl("-8$", colnames(sobj)))),
  rep("E18.5", sum(grepl("-9$", colnames(sobj))))
)

my_levels <- c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5", "E18.5")
sobj$stages <- factor(x = sobj$stages, levels = my_levels)

# sobj <- subset(sobj, subset = nFeature_RNA > 500 & percent_mt < 10) # no use, already filtered

sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$stages)
sobj

# run standard anlaysis workflow
sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj)

saveRDS(sobj, file = outRDS)
