#Packages to load
library(Mus.musculus)
library(magrittr)
library(org.Mm.eg.db)
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(data.table)
library(ComplexHeatmap)
library(slingshot)
library(scales)
library(scID)
library(umap)
library(tradeSeq)
library(BiocParallel)

set.seed(123)

B6NKT.data <- Read10X(data.dir = "/Volumes/Samsung_T5/CellRangerOutput/B6_1_total/filtered_feature_bc_matrix")
B6NKT <- CreateSeuratObject(counts = B6NKT.data, project = "B6_Thymus_NKT1")

B6_2NKT.data <- Read10X(data.dir = "/Volumes/Samsung_T5/CellRangerOutput/B6_2_total/filtered_feature_bc_matrix")
B6_2NKT <- CreateSeuratObject(counts = B6_2NKT.data, project = "B6_Thymus_NKT2")

Hivep31.exp.data <- Read10X(data.dir = "/Volumes/Samsung_T5/CellRangerOutput/Hivep3_1_total/filtered_feature_bc_matrix")
Hivep3_1 <- CreateSeuratObject(counts = Hivep31.exp.data, project = "Hivep3_Thymus_NKT1")
Hivep32.exp.data <- Read10X(data.dir = "/Volumes/Samsung_T5/CellRangerOutput/Hivep3_2_total/filtered_feature_bc_matrix")
Hivep3_2 <- CreateSeuratObject(counts = Hivep32.exp.data, project = "Hivep3_Thymus_NKT2")

Combined <- merge(B6NKT, y = c(B6_2NKT, Hivep3_1, Hivep3_2)
                add.cell.ids = c('B61', 'B62', 'Hivep31', 'Hivep32'), 
                project = "NKT")
Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern = "^mt-")
Combined <- subset(Combined, subset = nFeature_RNA >= 500 & nFeature_RNA < 4100 & nCount_RNA >=1000 & percent.mt < 6)

Combined <- NormalizeData(Combined)
Combined <- FindVariableFeatures(Combined)
Combined <- RunFastMNN(object.list = SplitObject(Combined, split.by = "orig.ident"), k = 20)
Combined <- RunUMAP(Combined, reduction = "mnn", dims = 1:30)
Combined <- FindNeighbors(Combined, reduction = "mnn", dims = 1:30)
Combined <- FindClusters(Combined, resolution = 0.8)

#Display UMAP
DimPlot(Combined, group.by = c("orig.ident", "seurat_clusters"), pt.size = 0.1, ncol = 2, label = TRUE)
DimPlot(Combined, group.by = c("seurat_clusters"), pt.size = 0.1, label = TRUE)
DimPlot(Combined, pt.size = 0.1, label = TRUE, split.by = "orig.ident", ncol = 4)

#Display differentially expressed genes on UMAP
#Display different genes by providing different list of genes to features parameter
FeaturePlot(Combined, features = c('Cd24a', 'Mki67', 'Zbtb16', 'Rorc', 'Tbx21'),
            cols = c('white', 'red'), pt.size = 0.7, ncol = 4)

#Create a heatmap of differentially expressed genes by cluster
diff.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- diff.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(pbmc, features = top5$gene) + NoLegend()

#Create a dotplot of genes highly expressed in different iNKT subsets
#Display different genes by providing different list of genes to markers.to.plot variable
markers.to.plot <- c("Cd24a", "Cd69", "Egr2", "Il2rb", "Ifng", "Tbx21", "Ccr6", "Rorc", "Il17rb", 
                     "Ccr7", "Il4", "Gata3", "Zbtb16", "Cd4", "Cd44")
DotPlot(Combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8)

#Create a violin plot for gene expression split up by sample
#Display different genes by providing different list of genes to features parameter
VlnPlot(Combined, features = c('Cd52', 'Drosha', 'Dusp1', 'Zbtb16', 'Plac8', 'Ccr9',
                               'Hdac7', 'Egr2', 'Egr1', 'Id3', 'Trac', 'Cd247', 'Trbc2', 'Trbc1'), pt.size = 0)


###################### PSEUDOTEMPORAL ANALYSIS USING SLINGSHOT #######################
#this helps match cluster colors from seurat UMAP to slingshot UMAP
cell_pal <- function(cell_vars, pal_fun, ...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
NKT.sce <- as.SingleCellExperiment(Combined)

#Perform slingshot analysis
NKT.sce <- slingshot(NKT.sce, reducedDim = 'UMAP', clusterLabels = 'seurat_clusters', start.clus = '0')
NKT.slo <- SlingshotDataSet(NKT.sce)

#Display trajectories on the UMAP
cell_colors_clust <- cell_pal(Idents(Combined), hue_pal())
plot(reducedDims(NKT.sce)$UMAP, col = cell_colors_clust, pch=16)
lines(SlingshotDataSet(NKT.sce), lwd=2, col='black')





