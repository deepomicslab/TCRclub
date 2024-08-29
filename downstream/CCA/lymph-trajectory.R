library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratWrappers)
#library(monocle3)
devtools::load_all("D:/TCRdataset/CCA_GSE201425/monocle_2.30.0/monocle")
library(patchwork)
library(umap)
library(anndata)
library(dplyr)
library(ggplot2)

datapath <- "D:/TCRdataset/CCA_GSE201425/myresult/integrate_filtered"
resultpath <- "D:/TCRdataset/CCA_GSE201425/myresult/integrate_filtered"
countpath <- "D:/TCRdataset/CCA_GSE201425/sample_filtered"
combined_filtered <- readRDS(file.path("D:/TCRdataset/CCA_GSE201425", "r_save", "LPcombined_filtered.rds"))
datadirs <- list.files(path=resultpath, full.names = FALSE)
datadirs <- datadirs[grep("LP", datadirs)]
lymph_node.list <-list()
for (dirname in datadirs){
  print(dirname)
  consensus <- read.table(file=file.path(resultpath, dirname, "consensus_result.csv"), header = TRUE, row.names="barcode", sep=",")
  countData <- read.table(file=file.path(countpath, sub("LP", "lymph_node_filtered", dirname), "rnacount.csv"), header = TRUE, row.names="barcode", sep=",")
  lymph_node <- CreateSeuratObject(counts = t(countData), project = sub("LP", "lymph_node_filtered", dirname), min.cells = 3, min.features = 200)
  lymph_node <- RenameCells(obj=lymph_node, add.cell.id = sub("LP", "lymph_node", dirname))
  #合并consensus result的PB_sort和obj@meta.data
  new_df <- data.frame(row.names=rownames(consensus),lymph_per=consensus$lymph_per, celltype=consensus$predicted_labels)
  lymph_node@meta.data <- merge(lymph_node@meta.data, new_df, by = "row.names", all.x = TRUE)
  row.names(lymph_node@meta.data) <- lymph_node@meta.data$Row.names
  lymph_node@meta.data$Row.names <- NULL
  
  lymph_node.list <- c(lymph_node.list, lymph_node)
  names(lymph_node.list)[length(lymph_node.list)] <- sub("LP", "lymph_node_filtered", dirname)
}
# normalize and identify variable features for each dataset independently
lymph_node.list <- lapply(X = lymph_node.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = lymph_node.list)
immune.anchors <- FindIntegrationAnchors(object.list = lymph_node.list, anchor.features = features)
lymph_nodeIntegrate <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(lymph_nodeIntegrate) <- "integrated"
lymph_nodeIntegrate <- ScaleData(lymph_nodeIntegrate, verbose = FALSE)
lymph_nodeIntegrate <- RunPCA(lymph_nodeIntegrate, npcs = 30, verbose = FALSE)
lymph_nodeIntegrate <- RunUMAP(lymph_nodeIntegrate, reduction = "pca", dims = 1:30, reduction.name = "umap")
lymph_nodeIntegrate <- FindNeighbors(lymph_nodeIntegrate,verbose = FALSE)
lymph_nodeIntegrate <- FindClusters(lymph_nodeIntegrate,verbose = FALSE, resolution = 0.5)
Idents(object = lymph_nodeIntegrate) <- "seurat_clusters"
p1 <- DimPlot(lymph_nodeIntegrate, reduction = "umap", label = TRUE, repel = TRUE)
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/Annotaionlymph_clusters.pdf", p1, width = 8, height = 6, dpi = 300)
DefaultAssay(lymph_nodeIntegrate) <- "RNA"
lymph_nodeIntegrate <- ScaleData(lymph_nodeIntegrate, verbose = FALSE)
features <- c("CD4", "CD8A", "CCR7", "IL7R", "TCF7","IL2RA","FOXP3",
              "PDCD1","CTLA4","LAG3","TIGIT","HAVCR2")
p1 <- VlnPlot(lymph_nodeIntegrate, features = features, pt.size=0) #小提琴图
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/Annotaionlymph_feauture.pdf", p1, width = 10, height = 10, dpi = 300)

features <- c("CD4", "CD8A", "CCR7", "IL7R", "TCF7","IL2RA","FOXP3",
              "PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","KLRD1","GZMK","GZMB","GNLY")
p1 <- DoHeatmap(lymph_nodeIntegrate, features = features) 
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/Annotaionlymph_heatmap.pdf", p1, width = 8, height = 6, dpi = 300)


#注释
new.cluster.ids <- c("Memory CD8 T cell", "Naive CD8 T cell", "Treg CD4", "Exhausted CD8 T cell", 
                     "NK T & Others","NK T & Others","NK T & Others","Exhausted CD8 T cell")
names(new.cluster.ids) <- levels(lymph_nodeIntegrate)
lymph_nodeIntegrate <- RenameIdents(lymph_nodeIntegrate, new.cluster.ids)
#在metadata中，添加Celltype信息
lymph_nodeIntegrate$celltype <- Idents(lymph_nodeIntegrate)
Idents(object = lymph_nodeIntegrate) <- "celltype"
p1 <- DimPlot(lymph_nodeIntegrate, reduction = "umap", label = TRUE, repel = TRUE)
p1
p1 <- DimPlot(lymph_nodeIntegrate, reduction = "umap", label = TRUE, repel = TRUE,
              cols = c('Memory CD8 T cell'='#F3627B', 'Naive CD8 T cell' = '#FF7F00','Treg CD4' = '#A65628',
                       "Exhausted CD8 T cell"="#33A02C", "NK T & Others"="#984EA3", "Naive CD4 T cell"="#CAB2D6", 
                       "NK T cell"="#6A3D9A","Effector CD8 T cell"="#78c6eb", "Helper CD4 T cell"="#d2993f"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/lymph_celltype.pdf", p1, width = 8, height = 5, dpi = 300)


lymph_nodeIntegrate <- subset(lymph_nodeIntegrate, subset = (celltype == c("Memory CD8 T cell","Exhausted CD8 T cell",
                                                               "NK T & Others","Effector CD8 T cell", "Naive CD8 T cell")))
lymph_nodeIntegrate <- subset(lymph_nodeIntegrate, subset = (celltype != "Treg CD4"))


FeaturePlot(lymph_nodeIntegrate, features = "lymph_per")
cellsfirst <- rownames(lymph_nodeIntegrate@meta.data[(lymph_nodeIntegrate@meta.data$lymph_per < 1),])
p1<-DimPlot(object = lymph_nodeIntegrate, cells.highlight = cellsfirst, cols.highlight = "#A61C5D", cols = "gray", order = TRUE)
cellsthird <- rownames(lymph_nodeIntegrate@meta.data[(lymph_nodeIntegrate@meta.data$lymph_per == 1),])
p3<-DimPlot(object = lymph_nodeIntegrate, cells.highlight = cellsthird, cols.highlight = "#1B64A4", cols = "grey", order = TRUE)
p3+p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/lymph_cellchange.pdf", p3+p1, width = 8, height = 5, dpi = 300)

annotationfile <- data.frame(barcode=rownames(lymph_nodeIntegrate@meta.data),lymph_nodeIntegrate@meta.data)
write.table(annotationfile, "D:/dataset/CCA_GSE201425/lymph_filtered_annotation.csv",  
            sep = ',', row.names = F, col.names = T, quote = F)
