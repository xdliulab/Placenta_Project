library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)

sobj <- readRDS("../data/sobj_x9_allcelltype_mnn.rds")
md <- sobj@meta.data
total <- md[ , c("orig.ident", "celltype")]
total$Barcode <- row.names(total)
total <- total[ , c("Barcode", "celltype")]
colnames(total) <- c("Barcode", "total_celltype")
precursor <- read.csv("../data/P-TGC-RNA_snn_res.1.4.csv")
colnames(precursor) <- c("Barcode", "precursor_celltype")
unique(precursor$precursor_celltype)
# # Assuming 'precursor' is your data frame
# precursor <- precursor %>%
#   filter(precursor_celltype == "DSC precursor")
# 
setdiff(precursor$Barcode, total$Barcode)

merged_data <- left_join(total, precursor, by = "Barcode")

# Replace values in total
merged_data <- merged_data %>%
  mutate(total_celltype = ifelse(!is.na(precursor_celltype), precursor_celltype, total_celltype))

# Selecting relevant columns
total <- merged_data %>% select(Barcode, total_celltype)
head(total)
unique(total$total_celltype)

all.equal(colnames(sobj), total$Barcode)
sobj$celltype <- total$total_celltype
Idents(sobj) <- "celltype"
DimPlot(sobj)
celltype <- unique(sobj$celltype)

color <- read.csv("../data/color_celltype_240125.csv", row.names = 1)
setdiff(color$order, celltype)
setdiff(celltype, color$order)

Idents(sobj) <- "celltype"
DimPlot(sobj, label = TRUE)

# Set the levels of the identity classes in your Seurat object
levels(sobj) <- color$order
color_vector <- setNames(color$Color, color$order)
p <- DimPlot(sobj, label = TRUE, cols = color_vector)
ggsave(filename = "../plots/UMAP_allcelltype_v3.pdf", plot = p, 
       device = "pdf", width = 14, height = 10)

saveRDS(sobj, file = "../data/sobj_x9_allcelltype_mnn_v3.rds")

