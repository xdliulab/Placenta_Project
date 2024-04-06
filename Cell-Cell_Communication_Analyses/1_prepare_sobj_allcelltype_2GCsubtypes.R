library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(loupeR)

job_id <- "627"

sobj <- readRDS("../../2024-03-15_620-x9-snRNA-allcelltype-annotation-v4/data/sobj_x9_allcelltype_mnn_v4.rds")
md <- sobj@meta.data
Idents(sobj) <- "celltype"
UMAPPlot(sobj, label = TRUE)

TB <- read.csv("../../2024-03-14_611-update-x9-TB-annotation-GC-subtype/output/TB_total_updated_v5.csv")
colnames(TB) <- c("Barcode", "celltype")

filtered_TB <- subset(TB, celltype %in% c("GC-Prl7b1", "GC-Aldh1a3"))

total <- data.frame(Barcode = row.names(md), celltype = md$celltype)
setdiff(filtered_TB$Barcode, total$Barcode)

merged_data <- left_join(total, filtered_TB, by = "Barcode")

# Replace values in total
merged_data <- merged_data %>%
  mutate(celltype.x = ifelse(!is.na(celltype.y), celltype.y, celltype.x))

total <- merged_data %>% select(Barcode, celltype.x)

all.equal(colnames(sobj), total$Barcode)

sobj$celltype <- total$celltype.x
Idents(sobj) <- "celltype"
UMAPPlot(sobj, label = TRUE)

clusters <- select_clusters(sobj)
projections <- select_projections(sobj)
count_mat <- LayerData(sobj, assay = "RNA", layer = "counts")
create_loupe(count_mat = count_mat, clusters = clusters, projections = projections, output_name = paste0("../output/loupe_snRNA_x9_allcelltype_2GCsubtypes_", job_id), force = TRUE)

saveRDS(sobj, file = paste0("../output/sobj_snRNA_x9_allcelltype_2GCsubtypes_", job_id, ".rds"))

sobj <- readRDS("../output/sobj_snRNA_x9_allcelltype_2GCsubtypes_627.rds")

celltype <- unique(sobj$celltype)

color <- read.csv("../data/color_celltype_240318.csv", row.names = 1)
setdiff(color$order, celltype)
setdiff(celltype, color$order)

# Set the levels of the identity classes in your Seurat object
DimPlot(sobj, label = TRUE)
levels(sobj) <- color$order
color_vector <- setNames(color$Color, color$order)
p <- DimPlot(sobj, label = TRUE, cols = color_vector)
ggsave(filename = "../plots/UMAP_allcelltype_627.pdf", plot = p, 
       device = "pdf", width = 14, height = 10)
