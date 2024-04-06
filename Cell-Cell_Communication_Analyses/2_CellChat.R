library(Seurat)
library(CellChat)
library(patchwork)
library(tidyverse)

sobj <- readRDS("../output/sobj_snRNA_x9_allcelltype_2GCsubtypes_627.rds")
# ref <- subset(sobj, subset = stage == input_stage)

cellchat <- createCellChat(object = sobj, group.by = "celltype", assay = "RNA")
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 15) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# execution.time = Sys.time() - ptm
# print(as.numeric(execution.time, units = "secs"))

options(future.globals.maxSize = 3 * 1024^3)  # Set limit to 3 GiB
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

saveRDS(df.net, file = "../data/df.net.rds")
saveRDS(cellchat, file = "../data/cellchat_snRNA_x9_allcelltype_2GCsubtypes.rds")
