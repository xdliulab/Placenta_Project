library(RColorBrewer)
library(SCENIC)
library(doParallel)
library(dplyr)
library(AUCell)
library(ComplexHeatmap)

scenicOptions <- readRDS("int/scenicOptions.Rds")

cellInfo <- readRDS("int/cellInfo.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists="null", verbose=FALSE)
cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
cellInfo <- data.frame(cellInfo)
colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")

binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

new_order <- c("JZP", "GC precursor", "GC-Aldh1a3", "GC-Prl7b1")
setdiff(new_order, colnames(regulonActivity_byCellType_Scaled))
setdiff(colnames(regulonActivity_byCellType_Scaled), new_order)
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[ , new_order]

pdf(file = "SCENIC_GC-Aldh1a3_reorder.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 48) # The height of the plot in inches

Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
        cluster_columns = FALSE)

dev.off()

row.names(regulonActivity_byCellType_Scaled)

filter_rows <- c( "Tfap2c (23g)", "Tcf12_extended (12g)", "Rora (15g)", "Rcor1 (10g)", "Nfyc_extended (14g)", "Nfib_extended (11g)", "Kdm2b (13g)", "Etv6 (11g)", "Ep300 (16g)", "Elf4_extended (30g)", "Cux1 (11g)", "Arid5b (11g)", "Anxa11 (57g)")

regulonActivity_byCellType_Scaled_filter <- regulonActivity_byCellType_Scaled[filter_rows, ]

pdf(file = "SCENIC-GC-Aldh1a3.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 6) # The height of the plot in inches

Heatmap(regulonActivity_byCellType_Scaled_filter, name="Regulon activity",
        cluster_columns = FALSE)

dev.off()

regulon <- readRDS("./int/2.6_regulons_asGeneSet.Rds")





