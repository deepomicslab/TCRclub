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
library(EnhancedVolcano)


#Integrate
resultpath <- "D:/TCRdataset/CCA_GSE201425/myresult/integrate_filtered"
countpath <- "D:/TCRdataset/CCA_GSE201425/sample_filtered"
combined_filtered <- readRDS(file.path("D:/TCRdataset/CCA_GSE201425", "r_save", "combined_filtered.rds"))
datadirs <- list.files(path=resultpath, full.names = FALSE)
datadirs <- datadirs[grep("primary", datadirs)]
primary.list <-list()
for (dirname in datadirs){
  print(dirname)
  consensus <- read.table(file=file.path(resultpath, dirname, "consensus_result.csv"), header = TRUE, row.names="barcode", sep=",")
  countData <- read.table(file=file.path(countpath, sub("PBMCprimary", "primary_focus_filtered", dirname), "rnacount.csv"), header = TRUE, row.names="barcode", sep=",")
  primary <- CreateSeuratObject(counts = t(countData), project = sub("PBMCprimary", "primary_focus_filtered", dirname), min.cells = 3, min.features = 200)
  primary <- RenameCells(obj=primary, add.cell.id = sub("PBMCprimary", "primary_focus", dirname))
  #合并consensus result的PB_sort和obj@meta.data
  new_df <- data.frame(row.names=rownames(consensus),PB_per=consensus$PB_per,celltype=consensus$predicted_labels)
  primary@meta.data <- merge(primary@meta.data, new_df, by = "row.names", all.x = TRUE)
  row.names(primary@meta.data) <- primary@meta.data$Row.names
  primary@meta.data$Row.names <- NULL
  
  primary.list <- c(primary.list, primary)
  names(primary.list)[length(primary.list)] <- sub("primary", "", dirname)
}
# normalize and identify variable features for each dataset independently
primary.list <- lapply(X = primary.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#primaryIntegrate <- merge(primary.list[[1]], primary.list[2:length(primary.list)], add.cell.ids = names(primary.list))
#primary.list <- SplitObject(primaryIntegrate, split.by = "orig.ident")
features <- SelectIntegrationFeatures(object.list = primary.list)
immune.anchors <- FindIntegrationAnchors(object.list = primary.list, anchor.features = features)
primaryIntegrate <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(primaryIntegrate) <- "integrated"
primaryIntegrate <- ScaleData(primaryIntegrate, verbose = FALSE)
primaryIntegrate <- RunPCA(primaryIntegrate, npcs = 30, verbose = FALSE)
primaryIntegrate <- RunUMAP(primaryIntegrate, reduction = "pca", dims = 1:30, reduction.name = "umap")
primaryIntegrate <- FindNeighbors(primaryIntegrate,verbose = FALSE)
primaryIntegrate <- FindClusters(primaryIntegrate,verbose = FALSE, resolution = 0.5)
Idents(object = primaryIntegrate) <- "seurat_clusters"
p1 <- DimPlot(primaryIntegrate, reduction = "umap", label = TRUE, repel = TRUE)
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/Annotaionprimary_clusters.pdf", p1, width = 8, height = 6, dpi = 300)
DefaultAssay(primaryIntegrate) <- "RNA"
features <- c("CD4", "CD8A", "CCR7", "IL7R", "TCF7","IL2RA","FOXP3",
              "PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","CCL5")
p1 <- VlnPlot(primaryIntegrate, features = features, pt.size=0) #小提琴图
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/Annotationprimary_feature.jpg", p1, width = 10, height = 10, dpi = 300)

primaryIntegrate <- ScaleData(primaryIntegrate, verbose = FALSE)
features <- c("CD4", "CD8A", "CCR7", "IL7R", "TCF7","IL2RA","FOXP3",
              "PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","KLRD1","GZMK","GZMB","GNLY","GZMH","IFITM1","CCL5")
p1 <- DoHeatmap(primaryIntegrate, features = features) 
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/Annotaionprimary_heatmap.jpg", p1, width = 8, height = 6, dpi = 300)

DefaultAssay(primaryIntegrate) <- "RNA"
genemarkers <- FindAllMarkers(primaryIntegrate, only.pos = TRUE)
genemarkers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() -> top10
DoHeatmap(primaryIntegrate, features = top10$gene) + NoLegend()
P1 <- DoHeatmap(primaryIntegrate, features=features, group.by="seurat_clusters")
P1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/Annotaionprimary_heatmap.pdf", p1, width = 8, height = 6, dpi = 300)

Idents(object = primaryIntegrate) <- "celltype"
p1 <- DimPlot(primaryIntegrate, reduction = "umap", label = TRUE, repel = TRUE)
p1



new.cluster.ids <- c("Effector CD8 T cell", "Exhausted CD8 T cell", "Memory CD8 T cell", "Treg CD4",
                     "Exhausted CD8 T cell", "Effector CD8 T cell", "Memory CD8 T cell", "Helper CD4 T cell",
                     "NK T cell") 
names(new.cluster.ids) <- levels(primaryIntegrate)
primaryIntegrate <- RenameIdents(primaryIntegrate, new.cluster.ids)
#在metadata中，添加Celltype信息
primaryIntegrate$celltype <- Idents(primaryIntegrate)
Idents(object = primaryIntegrate) <- "celltype"
p1 <- DimPlot(primaryIntegrate, reduction = "umap", label = TRUE, repel = TRUE,
              cols = c('Memory CD8 T cell'='#F3627B', 'Naive CD8 T cell' = '#FF7F00','Treg CD4' = '#A65628',
                       "Exhausted CD8 T cell"="#33A02C", "DN T & Others"="#984EA3", "Naive CD4 T cell"="#CAB2D6", 
                       "NK T cell"="#6A3D9A","Effector CD8 T cell"="#69A3DD", "Helper CD4 T cell"="#d2993f"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/primary_celltype2.pdf", p1, width = 8, height = 5, dpi = 300)

DimPlot(primaryIntegrate, group.by = c("seurat_clusters", "celltype"),reduction = "umap")


FeaturePlot(primaryIntegrate, features = "PB_per")
cellsfirst <- rownames(primaryIntegrate@meta.data[(primaryIntegrate@meta.data$PB_per == 0),])
p1<-DimPlot(object = primaryIntegrate, cells.highlight = cellsfirst, cols.highlight = "#484EAA", cols = "gray", order = TRUE)
cellsthird <- rownames(primaryIntegrate@meta.data[(primaryIntegrate@meta.data$PB_per > 0),])
p3<-DimPlot(object = primaryIntegrate, cells.highlight = cellsthird, cols.highlight = "#EC2835", cols = "grey", order = TRUE)
p3+p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/primary_cellchange.pdf", p3+p1, width = 8, height = 5, dpi = 300)

#Effector
effecttumor <- subset(primaryIntegrate, subset = celltype == c("Effector CD8 T cell"))
DefaultAssay(effecttumor) <- "RNA"
Idents(object = effecttumor) <- "celltype"
DimPlot(effecttumor, reduction = "umap", label = TRUE, repel = TRUE)

cellsfirst <- rownames(effecttumor@meta.data[(effecttumor@meta.data$PB_per ==0)&(effecttumor@meta.data$celltype == "Effector CD8 T cell"),])
p1<-DimPlot(object = effecttumor, cells.highlight = cellsfirst, cols.highlight = "red", cols = "gray", order = TRUE)
cellssecond <- rownames(effecttumor@meta.data[(effecttumor@meta.data$PB_per > 0)&(effecttumor@meta.data$celltype == "Effector CD8 T cell"),])
p2<-DimPlot(object = effecttumor, cells.highlight = cellssecond, cols.highlight = "red", cols = "gray", order = TRUE)
p1+p2
effecttumor@meta.data$stage <- "None"
effecttumor@meta.data[cellsfirst,]$stage <- "low"
effecttumor@meta.data[cellssecond,]$stage <- "high"
Idents(object = effecttumor) <- "stage"
deg.cluster <- FindMarkers(effecttumor, ident.1="low", ident.2="high")
deg.cluster$p_val_adj <- p.adjust(p=deg.cluster[,1],method="BH", n=nrow(effecttumor))
deg.cluster <- deg.cluster[deg.cluster$p_val_adj<0.05,]
p1<-EnhancedVolcano(deg.cluster,
                    lab=rownames(deg.cluster), x = 'avg_log2FC', y = 'p_val_adj',
                    gridlines.major=FALSE, gridlines.minor=FALSE,pCutoff = 0.05,FCcutoff = 0.5)
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/DE/EffectorVolcano2.pdf", p1, width = 6, height = 5, dpi = 300)

#Effector
effecttumor <- subset(primaryIntegrate, subset = celltype == c("Effector CD8 T cell", "Exhausted CD8 T cell"))
DefaultAssay(effecttumor) <- "RNA"
Idents(object = effecttumor) <- "celltype"
DimPlot(effecttumor, reduction = "umap", label = TRUE, repel = TRUE)

cellsfirst <- rownames(effecttumor@meta.data[(effecttumor@meta.data$PB_per ==0)&(effecttumor@meta.data$celltype == "Effector CD8 T cell"),])
p1<-DimPlot(object = effecttumor, cells.highlight = cellsfirst, cols.highlight = "red", cols = "gray", order = TRUE)
cellssecond <- rownames(effecttumor@meta.data[(effecttumor@meta.data$PB_per > 0)&(effecttumor@meta.data$celltype == "Effector CD8 T cell"),])
p2<-DimPlot(object = effecttumor, cells.highlight = cellssecond, cols.highlight = "red", cols = "gray", order = TRUE)
p1+p2
effecttumor@meta.data$stage <- "Exhausted CD8 T cell"
effecttumor@meta.data[cellsfirst,]$stage <- "low"
effecttumor@meta.data[cellssecond,]$stage <- "high"
Idents(object = effecttumor) <- "stage"
deg.cluster <- FindMarkers(effecttumor, ident.1="low", ident.2="high")
deg.cluster$p_val_adj <- p.adjust(p=deg.cluster[,1],method="BH", n=nrow(effecttumor))
deg.cluster <- deg.cluster[deg.cluster$p_val_adj<0.05,]
express_genes <- rownames(deg.cluster)

#Plot Trajectory
data <- as(as.matrix(effecttumor@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = effecttumor@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       #lowerDetectionLimit = 0.5,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150
HSMM <- setOrderingFilter(HSMM, express_genes)
HSMM <- reduceDimension(HSMM, norm_method="none", max_components = 3, method = "DDRTree") # 可调参数，选择降维方法和主成分
HSMM <- orderCells(HSMM)
p1<-plot_cell_trajectory(HSMM, color_by = "stage")+ scale_color_manual(values = c("gray","red", "blue"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Effector_HighLow2.pdf", p1, width = 6, height = 5, dpi = 300)
p1 <- plot_cell_trajectory(HSMM, color_by = "State")
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Effector_state.pdf", p1, width = 6, height = 5, dpi = 300)

table(HSMM@phenoData@data[HSMM@phenoData@data$State == 1,]$stage)
write.csv(HSMM@phenoData@data, "Effector_state.csv")

p1<-plot_cell_trajectory(HSMM, color_by = "stage") + facet_wrap(~State, nrow = 1)+ scale_color_manual(values = c("gray","red", "blue"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/EffectorHighLow_resolution.pdf", p1, width = 6, height = 5, dpi = 300)
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HSMM_1 <- orderCells(HSMM, root_state = GM_state(HSMM))
p1<-plot_cell_trajectory(HSMM_1, color_by = "Pseudotime")
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Effectortrajectory2.pdf", p1, width = 6, height = 5, dpi = 300)

#Memory
memorytumor <- subset(primaryIntegrate, subset = celltype == c("Memory CD8 T cell"))
DefaultAssay(memorytumor) <- "RNA"
Idents(object = memorytumor) <- "celltype"
DimPlot(memorytumor, reduction = "umap", label = TRUE, repel = TRUE)

cellsfirst <- rownames(memorytumor@meta.data[(memorytumor@meta.data$PB_per == 0)&(memorytumor@meta.data$celltype == "Memory CD8 T cell"),])
p1<-DimPlot(object = memorytumor, cells.highlight = cellsfirst, cols.highlight = "red", cols = "gray", order = TRUE)
cellssecond <- rownames(memorytumor@meta.data[(memorytumor@meta.data$PB_per >0)&(memorytumor@meta.data$celltype == "Memory CD8 T cell"),])
p2<-DimPlot(object = memorytumor, cells.highlight = cellssecond, cols.highlight = "red", cols = "gray", order = TRUE)
p2+p1
memorytumor@meta.data$stage <- "None"
memorytumor@meta.data[cellsfirst,]$stage <- "low"
memorytumor@meta.data[cellssecond,]$stage <- "high"
Idents(object = memorytumor) <- "stage"
DefaultAssay(memorytumor) <- "RNA"
deg.cluster <- FindMarkers(memorytumor, ident.1="low", ident.2="high")
deg.cluster$p_val_adj <- p.adjust(p=deg.cluster[,1],method="BH", n=nrow(memorytumor))
deg.cluster <- deg.cluster[deg.cluster$p_val_adj<0.05,]
p1<-EnhancedVolcano(deg.cluster,
                    lab=rownames(deg.cluster), x = 'avg_log2FC', y = 'p_val_adj',
                    gridlines.major=FALSE, gridlines.minor=FALSE, pCutoff = 0.05,FCcutoff = 0.5)
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/DE/MemoryVolcano2.pdf", p1, width = 6, height = 5, dpi = 300)


memorytumor <- subset(primaryIntegrate, subset = celltype == c("Memory CD8 T cell", "Exhausted CD8 T cell"))
DefaultAssay(memorytumor) <- "RNA"
Idents(object = memorytumor) <- "celltype"
DimPlot(memorytumor, reduction = "umap", label = TRUE, repel = TRUE)
cellsfirst <- rownames(memorytumor@meta.data[(memorytumor@meta.data$PB_per == 0)&(memorytumor@meta.data$celltype == "Memory CD8 T cell"),])
p1<-DimPlot(object = memorytumor, cells.highlight = cellsfirst, cols.highlight = "red", cols = "gray", order = TRUE)
cellssecond <- rownames(memorytumor@meta.data[(memorytumor@meta.data$PB_per > 0)&(memorytumor@meta.data$celltype == "Memory CD8 T cell"),])
p2<-DimPlot(object = memorytumor, cells.highlight = cellssecond, cols.highlight = "red", cols = "gray", order = TRUE)
p2+p1
memorytumor@meta.data$stage <- "Exhausted CD8 T cell"
memorytumor@meta.data[cellsfirst,]$stage <- "low"
memorytumor@meta.data[cellssecond,]$stage <- "high"
Idents(object = memorytumor) <- "stage"
DefaultAssay(memorytumor) <- "RNA"
deg.cluster <- FindMarkers(memorytumor, ident.1="low", ident.2="high")
deg.cluster$p_val_adj <- p.adjust(p=deg.cluster[,1],method="BH", n=nrow(memorytumor))
deg.cluster <- deg.cluster[deg.cluster$p_val_adj<0.05,]
express_genes <- rownames(deg.cluster)

##因express_genes太少 使用seurat选择的高变基因
memorytumor <- FindVariableFeatures(memorytumor, selection.method = "vst", nfeatures = 2000)
express_genes <- VariableFeatures(memorytumor)

data <- as(as.matrix(memorytumor@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = memorytumor@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       #lowerDetectionLimit = 0.5,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150
HSMM <- setOrderingFilter(HSMM, express_genes)
HSMM <- reduceDimension(HSMM, norm_method="none", max_components = 2, method = "DDRTree",auto_param_selection = F) # 可调参数，选择降维方法和主成分
HSMM <- orderCells(HSMM)
p1 <- plot_cell_trajectory(HSMM, color_by = "stage")+ scale_color_manual(values = c("gray","red", "blue"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Memory_HighLow2.pdf", p1, width = 6, height = 5, dpi = 300)
p1 <- plot_cell_trajectory(HSMM, color_by = "State")
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Memory_state.pdf", p1, width = 6, height = 5, dpi = 300)
write.csv(HSMM@phenoData@data, "Memory_state.csv")

table(HSMM@phenoData@data[HSMM@phenoData@data$State == 6,]$stage)

HSMM <- orderCells(HSMM, root_state = 6)
p1 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Memoryrajectory2.pdf", p1, width = 6, height = 5, dpi = 300)

p1<-plot_cell_trajectory(HSMM, color_by = "stage") + facet_wrap(~State, nrow = 1)+ scale_color_manual(values = c("gray","red", "blue"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/MemoryHighLow_resolution.pdf", p1, width = 6, height = 5, dpi = 300)


#Exhausted
exhaustedtumor <- subset(primaryIntegrate, subset = celltype == c("Exhausted CD8 T cell"))
DefaultAssay(exhaustedtumor) <- "RNA"
Idents(object = exhaustedtumor) <- "celltype"
DimPlot(exhaustedtumor, reduction = "umap", label = TRUE, repel = TRUE)

cellsfirst <- rownames(exhaustedtumor@meta.data[(exhaustedtumor@meta.data$PB_per == 0)&(exhaustedtumor@meta.data$celltype == "Exhausted CD8 T cell"),])
p1<-DimPlot(object = exhaustedtumor, cells.highlight = cellsfirst, cols.highlight = "red", cols = "gray", order = TRUE)
cellssecond <- rownames(exhaustedtumor@meta.data[(exhaustedtumor@meta.data$PB_per > 0)&(exhaustedtumor@meta.data$celltype == "Exhausted CD8 T cell"),])
p2<-DimPlot(object = exhaustedtumor, cells.highlight = cellssecond, cols.highlight = "red", cols = "gray", order = TRUE)
p2+p1
exhaustedtumor@meta.data$stage <- "None"
exhaustedtumor@meta.data[cellsfirst,]$stage <- "low"
exhaustedtumor@meta.data[cellssecond,]$stage <- "high"
Idents(object = exhaustedtumor) <- "stage"
DefaultAssay(exhaustedtumor) <- "RNA"
deg.cluster <- FindMarkers(exhaustedtumor, ident.1="low", ident.2="high")
deg.cluster$p_val_adj <- p.adjust(p=deg.cluster[,1],method="BH", n=nrow(exhaustedtumor))
deg.cluster <- deg.cluster[deg.cluster$p_val_adj<0.05,]
p1<-EnhancedVolcano(deg.cluster,
                    lab=rownames(deg.cluster), x = 'avg_log2FC', y = 'p_val_adj',
                    gridlines.major=FALSE, gridlines.minor=FALSE, pCutoff = 0.05,FCcutoff = 0.5)
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/DE/ExhaustedVolcano2.pdf", p1, width = 6, height = 5, dpi = 300)
express_genes <- rownames(subset(deg.cluster,p_val_adj<0.05))

##1使用seurat选择的高变基因⚠️
exhaustedtumor <- FindVariableFeatures(exhaustedtumor, selection.method = "vst", nfeatures = 2000)
express_genes <- VariableFeatures(exhaustedtumor)

data <- as(as.matrix(exhaustedtumor@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = exhaustedtumor@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       #lowerDetectionLimit = 0.5,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150
HSMM <- setOrderingFilter(HSMM, express_genes)
HSMM <- reduceDimension(HSMM, norm_method="none", max_components = 2, method = "DDRTree",auto_param_selection = F) # 可调参数，选择降维方法和主成分
HSMM <- orderCells(HSMM)
p1 <- plot_cell_trajectory(HSMM, color_by = "stage")+ scale_color_manual(values = c("red", "blue"))
p1
plot_cell_trajectory(HSMM, color_by = "State")
table(HSMM@phenoData@data[HSMM@phenoData@data$State == 3,]$stage)
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Exhausted_HighLow2.pdf", p1, width = 6, height = 5, dpi = 300)
p1<-plot_cell_trajectory(HSMM, color_by = "stage") + facet_wrap(~State, nrow = 1)+ scale_color_manual(values = c("gray","red", "blue"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/ExhaustedHighLow_resolution.pdf", p1, width = 6, height = 5, dpi = 300)
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HSMM_1 <- orderCells(HSMM, root_state = GM_state(HSMM))
p1<-plot_cell_trajectory(HSMM_1, color_by = "Pseudotime")
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/Exhaustedrajectory2.pdf", p1, width = 6, height = 5, dpi = 300)

#CD8tumor
CD8tumor <- subset(primaryIntegrate, subset = (celltype != c("Treg CD4")))
CD8tumor <- subset(CD8tumor, subset = (celltype != c("Helper CD4 T cell")))
Idents(object = CD8tumor) <- "celltype"

cellsfirst <- rownames(CD8tumor@meta.data[(CD8tumor@meta.data$PB_per <= 0.1)&(CD8tumor@meta.data$celltype == "Effector CD8 T cell"),])
p1<-DimPlot(object = CD8tumor, cells.highlight = cellsfirst, cols.highlight = "red", cols = "gray", order = TRUE)
cellssecond <- rownames(CD8tumor@meta.data[(CD8tumor@meta.data$PB_per > 0.1)&(CD8tumor@meta.data$celltype == "Effector CD8 T cell"),])
p2<-DimPlot(object = CD8tumor, cells.highlight = cellssecond, cols.highlight = "red", cols = "gray", order = TRUE)

cellsfirst <- rownames(CD8tumor@meta.data[(CD8tumor@meta.data$PB_per <= 0.1)&(CD8tumor@meta.data$celltype == "Memory CD8 T cell"),])
p1<-DimPlot(object = CD8tumor, cells.highlight = cellsfirst, cols.highlight = "red", cols = "gray", order = TRUE)
cellssecond <- rownames(CD8tumor@meta.data[(CD8tumor@meta.data$PB_per > 0.1)&(CD8tumor@meta.data$celltype == "Memory CD8 T cell"),])
p2<-DimPlot(object = CD8tumor, cells.highlight = cellssecond, cols.highlight = "red", cols = "gray", order = TRUE)

CD8tumor@meta.data$stage <- "Other Cells"
CD8tumor@meta.data[cellsfirst,]$stage <- "low"
CD8tumor@meta.data[cellssecond,]$stage <- "high"
Idents(object = CD8tumor) <- "stage"
DefaultAssay(CD8tumor) <- "RNA"
##1使用seurat选择的高变基因⚠️
express_genes <- VariableFeatures(CD8tumor)
data <- as(as.matrix(CD8tumor@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = CD8tumor@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       #lowerDetectionLimit = 0.5,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150
HSMM <- setOrderingFilter(HSMM, express_genes)
HSMM <- reduceDimension(HSMM, norm_method="none", max_components = 3, method = "DDRTree") # 可调参数，选择降维方法和主成分
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "stage")+ scale_color_manual(values = c("red", "blue","gray"))
plot_cell_trajectory(HSMM, color_by = "stage") + facet_wrap(~State, nrow = 1)
p1<-plot_cell_trajectory(HSMM, color_by = "celltype") + scale_color_manual(values=c("#69A3DD","#33A02C",'#F3627B',"#6A3D9A"))
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/CD8_trajectory.pdf", p1, width = 6, height = 5, dpi = 300)

plot_cell_trajectory(HSMM, color_by = "celltype")
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HSMM_1 <- orderCells(HSMM, root_state = GM_state(HSMM))
p1 <- plot_cell_trajectory(HSMM_1, color_by = "Pseudotime")
p1
ggsave("D:/论文/TCRclub Manuscript/figures/CCA/trajectory/CD8_pseudotime.pdf", p1, width = 6, height = 5, dpi = 300)


#导出注释文件
annotationfile <- data.frame(barcode=rownames(primaryIntegrate@meta.data),primaryIntegrate@meta.data)
write.table(annotationfile, "D:/dataset/CCA_GSE201425/primary_filtered_annotation.csv",  
            sep = ',', row.names = F, col.names = T, quote = F)


