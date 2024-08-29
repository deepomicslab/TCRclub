library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratWrappers)
#library(monocle3)
devtools::load_all("D:/dataset/CCA_GSE201425/monocle_2.30.0/monocle")
library(patchwork)
library(umap)
library(anndata)
library(dplyr)
library(ggplot2)

datapath <- "D:/dataset/CCA_GSE201425/myresult/integrate_filtered"


#Integrate版本
resultpath <- "D:/dataset/CCA_GSE201425/myresult/integrate_filtered"
countpath <- "D:/dataset/CCA_GSE201425/sample_filtered"
combined_filtered <- readRDS(file.path("D:/dataset/CCA_GSE201425", "r_save", "combined_filtered.rds"))
datadirs <- list.files(path=resultpath, full.names = FALSE)
datadirs <- datadirs[grep("PBMC", datadirs)]
PBMC.list <-list()
for (dirname in datadirs){
  print(dirname)
  consensus <- read.table(file=file.path(resultpath, dirname, "consensus_result.csv"), header = TRUE, row.names="barcode", sep=",")
  countData <- read.table(file=file.path(countpath, sub("primary", "_filtered", dirname), "rnacount.csv"), header = TRUE, row.names="barcode", sep=",")
  PBMC <- CreateSeuratObject(counts = t(countData), project = sub("primary", "_filtered", dirname), min.cells = 3, min.features = 200)
  PBMC <- RenameCells(obj=PBMC, add.cell.id = sub("primary", "", dirname))
  #合并consensus result的PB_sort和obj@meta.data
  new_df <- data.frame(row.names=rownames(consensus),PB_per=consensus$PB_per,
                       celltype=consensus$predicted_labels)
  PBMC@meta.data <- merge(PBMC@meta.data, new_df, by = "row.names", all.x = TRUE)
  row.names(PBMC@meta.data) <- PBMC@meta.data$Row.names
  PBMC@meta.data$Row.names <- NULL
  
  PBMC.list <- c(PBMC.list, PBMC)
  names(PBMC.list)[length(PBMC.list)] <- sub("primary", "", dirname)
}
# normalize and identify variable features for each dataset independently
PBMC.list <- lapply(X = PBMC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#PBMCIntegrate <- merge(PBMC.list[[1]], PBMC.list[2:length(PBMC.list)], add.cell.ids = names(PBMC.list))
#PBMC.list <- SplitObject(PBMCIntegrate, split.by = "orig.ident")
features <- SelectIntegrationFeatures(object.list = PBMC.list)
immune.anchors <- FindIntegrationAnchors(object.list = PBMC.list, anchor.features = features)
PBMCIntegrate <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(PBMCIntegrate) <- "integrated"
PBMCIntegrate <- ScaleData(PBMCIntegrate, verbose = FALSE)
PBMCIntegrate <- RunPCA(PBMCIntegrate, npcs = 30, verbose = FALSE)
PBMCIntegrate <- RunUMAP(PBMCIntegrate, reduction = "pca", dims = 1:30, reduction.name = "umap")
PBMCIntegrate <- FindNeighbors(PBMCIntegrate,verbose = FALSE)
PBMCIntegrate <- FindClusters(PBMCIntegrate,verbose = FALSE, resolution = 0.5)
Idents(object = PBMCIntegrate) <- "seurat_clusters"
p1 <- DimPlot(PBMCIntegrate, reduction = "umap", label = TRUE, repel = TRUE)
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/AnnotaionPBMC_clusters.pdf", p1, width = 8, height = 6, dpi = 300)
DefaultAssay(PBMCIntegrate) <- "RNA"
features <- c("CD4", "CD8A", "CCR7", "IL7R", "TCF7","IL2RA","FOXP3",
              "PDCD1","CTLA4","LAG3","TIGIT","HAVCR2")
p1 <- VlnPlot(PBMCIntegrate, features = features, pt.size=0) #小提琴图
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/AnnotaionPBMC_feauture.pdf", p1, width = 10, height = 10, dpi = 300)
PBMCIntegrate <- ScaleData(PBMCIntegrate, verbose = FALSE)
features <- c("CD4", "CD8A", "CCR7", "IL7R", "TCF7","IL2RA","FOXP3",
              "PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","KLRD1","GZMK","GZMB","GNLY")
p1 <- DoHeatmap(PBMCIntegrate, features = features) 
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/AnnotaionPBMC_heatmap.pdf", p1, width = 8, height = 6, dpi = 300)


Idents(object = PBMCIntegrate) <- "celltype"
p1 <- DimPlot(PBMCIntegrate, reduction = "umap", label = TRUE, repel = TRUE)
p1


#注释
new.cluster.ids <- c("Effector CD8 T cell", "Helper CD4 T cell", "Memory CD8 T cell", 
                     "Naive CD4 T cell", "Exhausted CD8 T cell", "NK T & Others","NK T & Others","NK T & Others")
names(new.cluster.ids) <- levels(PBMCIntegrate)
PBMCIntegrate <- RenameIdents(PBMCIntegrate, new.cluster.ids)
#在metadata中，添加Celltype信息
PBMCIntegrate$celltype <- Idents(PBMCIntegrate)
Idents(object = PBMCIntegrate) <- "celltype"
p1 <- DimPlot(PBMCIntegrate, reduction = "umap", label = TRUE, repel = TRUE,
              cols = c('Memory CD8 T cell'='#F3627B', 'Naive CD8 T cell' = '#FF7F00','Treg CD4' = '#A65628',
                       "Exhausted CD8 T cell"="#33A02C", "NK T & Others"="#984EA3", "Naive CD4 T cell"="#CAB2D6", 
                       "NK T cell"="#6A3D9A","Effector CD8 T cell"="#78c6eb", "Helper CD4 T cell"="#d2993f"))
p1
ggsave("D:/论文/TCRclub Manuscript/f>= 0.9igures/CCA/PBMC_celltype.pdf", p1, width = 8, height = 5, dpi = 300)

DimPlot(PBMCIntegrate, group.by = c("seurat_clusters", "celltype"),reduction = "umap")

PBMCIntegrate <- subset(PBMCIntegrate, subset = (celltype == c("Memory CD8 T cell","Exhausted CD8 T cell",
                                                               "NK T & Others","Effector CD8 T cell", "Naive CD8 T cell")))
PBMCIntegrate <- subset(PBMCIntegrate, subset = (celltype != "Naive CD4 T cell"))


FeaturePlot(PBMCIntegrate, features = "PB_per")
cellsfirst <- rownames(PBMCIntegrate@meta.data[(PBMCIntegrate@meta.data$PB_per <1),])
p1<-DimPlot(object = PBMCIntegrate, cells.highlight = cellsfirst, sizes.highlight=0.5, cols.highlight = "#A61C5D", cols = "gray", order = TRUE)
cellsthird <- rownames(PBMCIntegrate@meta.data[(PBMCIntegrate@meta.data$PB_per ==1),])
p3<-DimPlot(object = PBMCIntegrate, cells.highlight = cellsthird, sizes.highlight=0.5,cols.highlight = "#1B64A4", cols = "grey", order = TRUE)
p3+p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/PBMC_cellchange.pdf", p3+p1, width = 8, height = 5, dpi = 300)


#导出注释文件
annotationfile <- data.frame(barcode=rownames(PBMCIntegrate@meta.data),PBMCIntegrate@meta.data)
write.table(annotationfile, "D:/dataset/CCA_GSE201425/PBMC_filtered_annotation.csv",  
            sep = ',', row.names = F, col.names = T, quote = F)
