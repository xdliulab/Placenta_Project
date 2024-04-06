library(CellChat)
library(patchwork)
library(openxlsx)
library(reshape2)

dataLR <- readRDS("../data/df.net.rds")
write.xlsx(dataLR, file = "../output/CellChat_snRNA_x9_allcelltype_2GCsubtypes_integrate.xlsx")


###########
# Stop here
cellchat.HET <- readRDS("../data/cellchat_HET.rds")
dataLR_HET <- subsetCommunication(cellchat.HET)
write.xlsx(dataLR_HET, file = "../output/CellChat_HET.xlsx")

cellchat.KO <- readRDS("../data/cellchat_KO.rds")
dataLR_KO <- subsetCommunication(cellchat.KO)
write.xlsx(dataLR_KO, file = "../output/CellChat_KO.xlsx")

dataLR_WT$stage <- "WT"
dataLR_HET$stage <- "HET"
dataLR_KO$stage <- "KO"

dataLR <- rbind(dataLR_WT, dataLR_HET, dataLR_KO)
dataLR$Var1 <- dataLR$interaction_name_2
dataLR$Var2 <- paste0(dataLR$source, " - ", dataLR$target)
dataLR$significant <- dataLR$pval

data <- dataLR[ , c("Var1", "Var2", "significant", "stage")]
data$Var2_stage <- paste0(data$Var2, " (", data$stage, ")")
wide_data <- dcast(data, Var1 ~ Var2_stage, value.var = "significant", fill = "no")
column_order <- order(names(wide_data)[-1])
wide_data <- wide_data[, c(1, column_order + 1)]
write.csv(wide_data, file = "../data/heatmap_data.csv")
write.xlsx(wide_data, file = "../output/CellChat_heatmap_data.xlsx")


object.list <- list(KO = cellchat.KO, WT = cellchat.WT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
saveRDS(object.list, file = "../data/cellchat_object.list_KO_WT.rds")
saveRDS(cellchat, file = "../data/cellchat_merged_KO_WT.rds")



gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave("../plots/compareInteractions_KO_vs_WT.pdf", gg1, width = 7, height = 7, units = "in")

netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

netVisual_heatmap(cellchat)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

dataLR <- subsetCommunication(cellchat)


