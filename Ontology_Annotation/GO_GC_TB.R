library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stats)
library(data.table)
library(dplyr)
library(openxlsx)

all_markers <- readRDS("../data/DEG_TB_GC_subtypes.rds")

# Loop through each unique cluster in the data
unique_clusters <- unique(all_markers$cluster)

for (cluster_name in unique_clusters) {
  DEG_data <- all_markers %>% filter(cluster == cluster_name)
  
  # Gene name to GeneID conversion
  gene_df <- bitr(DEG_data$gene, fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Mm.eg.db)
  
  colnames(gene_df)[1] <- "gene"
  DEG_data1 <- left_join(gene_df, DEG_data)
  
  # GO enrichment
  GO_all <- enrichGO(gene = DEG_data1$ENTREZID,
                     keyType = "ENTREZID",
                     OrgDb = org.Mm.eg.db,
                     ont = "ALL",
                     pvalueCutoff = 0.01,
                     pAdjustMethod = "fdr",
                     minGSSize = 10,
                     maxGSSize = 500,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  
  GO_result <- data.frame(GO_all)
  
  # Write GO_result to Excel
  write.xlsx(GO_result, paste0("../output/GO_", cluster_name, ".xlsx"))
  
  # # Top 12 significant GO terms for plotting
  # go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 12, wt = -p.adjust)
  # 
  # # Create ggplot
  # p <- ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  #   geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  #   theme_minimal() +
  #   labs(x="GO Term", y="Gene_Number", title="Top 12 Enriched GO Terms") +
  #   facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y') +
  #   coord_flip() +
  #   scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen")) +
  #   theme_bw()
  # 
  # # Save ggplot to PDF
  # pdf_file <- paste0("../plots/", cluster_name, "_GO_plot.pdf")
  # pdf(pdf_file)
  # print(p)
  # dev.off()
}
