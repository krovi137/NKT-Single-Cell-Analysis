#Load required packages
library(tidyverse)
library(Seurat)
library(Matrix)
library(data.table)
library(scales)
library(viridis)
library(org.Mm.eg.db)
library(Mus.musculus)
library(msigdbr)
library(ComplexHeatmap)
library(EnrichmentBrowser)
library(clusterProfiler)
library(ggpmisc)
library(corrgram)
library(pheatmap)
library(RColorBrewer)
library(reshape)

#Load the required Seurat object and split into NKT and MAIT
B6.NKT.MAIT <- readRDS('/Volumes/Samsung_T5/NKT_Analysis/NKT_MAIT_Analysis/B6NKT_MAIT.rds')
Idents(B6.NKT.MAIT) <- 'orig.ident'
NKT <- subset(B6.NKT.MAIT, idents = c('B6NKT', 'B6_2NKT'))
Idents(NKT) <- 'New_Clusters'
MAIT <- subset(B6.NKT.MAIT, idents = 'MAIT')
Idents(MAIT) <- 'New_Clusters'

#Get average expression for each gene in each cluster for both NKT and MAIT samples
NKT.average.cluster <- AverageExpression(NKT, use.scale = TRUE, verbose = TRUE)[[1]]
colnames(NKT.average.cluster) <- paste0('NKT ', colnames(NKT.average.cluster))
MAIT.average.cluster <- AverageExpression(MAIT, use.scale = TRUE, verbose = TRUE)[[1]]
colnames(MAIT.average.cluster) <- paste0('MAIT ', colnames(MAIT.average.cluster))
Total.average <- as.matrix(cbind(NKT.average.cluster, MAIT.average.cluster))
Total.average <- Total.average[rowSums(Total.average) > 0, ]

#Determine Euclidean distances between the different clusters belonging to NKT and MAIT samples
euclid.dists <- dist(t(Total.average))
euclidMat <- as.matrix(euclid.dists)
euclidMat[upper.tri(euclidMat)] <- NA

euclidMat <- as.data.frame(euclidMat) %>%
  rownames_to_column(var = 'Samples')

#Reshape Euclidean matrix
euclidMat.reshape <- na.omit(melt(euclidMat, id=c('Samples'), variable_name = 'Comparators')) %>%
  dplyr::rename(Distance = value)
euclidMat.reshape$Samples <- factor(euclidMat.reshape$Samples, levels = levels(euclidMat.reshape$Comparators))
euclidMat.reshape$Comparators <- factor(euclidMat.reshape$Comparators, levels = rev(levels(euclidMat.reshape$Comparators)))

#Plot Euclidean matrix heatmap
a <- ggplot(euclidMat.reshape, aes(Samples, Comparators)) +
  theme_bw() +
  geom_tile(aes(fill = Distance), color='white') +
  coord_equal() +
  scale_fill_binned(breaks = c(20, 40, 60)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=270),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(position = "top")

print(a, vp = viewport(width = unit(8, "in"),
                       height = unit(8, "in"), angle = 45))