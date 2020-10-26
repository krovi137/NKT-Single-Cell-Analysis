#Packages to load
library(Mus.musculus)
library(magrittr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(EnrichmentBrowser)
library(biomaRt)
library(data.table)
library(ComplexHeatmap)
library(msigdbr)
library(scales)

#Function to identify enriched GO terms for each cluster
getEnrichedGO <- function(seuratObject, sig.db) {
  cluster.idents <- as.integer(levels(Idents(seuratObject)))
  gse.list <- lapply(X = cluster.idents, FUN = function(x) {
    print(paste0('Currently finding markers for cluster: ', x))
    diff.genes <- FindConservedMarkers(seuratObject, ident.1 = x, test.use = 'poisson', grouping.var = 'orig.ident',
                                       logfc.threshold = 0.25, min.pct = 0.25, only.pos = TRUE, verbose = TRUE,
                                       min.cells.group = 15)
    diff.genes$GeneID <- AnnotationDbi::select(Mus.musculus, keys=rownames(diff.genes), columns=c("ENTREZID"), keytype="SYMBOL")$ENTREZID
    diff.genes <- diff.genes[diff.genes$Hivep3NKT_p_val_adj < 0.01, ]
    gene.list <- diff.genes$GeneID
    
    em <- enricher(gene.list, TERM2GENE = sig.db, pAdjustMethod = 'BH',
                   minGSSize = 25, maxGSSize = 800, qvalueCutoff = 0.05)
    em.df <- em@result
    em.df <- em.df[em.df$qvalue <= 0.05, ]
    em.df
  })
  pathways <- rbindlist(gse.list, fill = TRUE, idcol = 'Cluster')
  pathways
}

#this is the database for gene sets that contain genes annotated by the same ontology term
#this will be used in the GO term enrichment analysis
C5_df <- msigdbr(species = "Mus musculus", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene) %>%
  mutate(gene = as.character(entrez_gene)) %>%
  dplyr::select(c(ont = gs_name, gene)) %>%
  as.data.frame()

#B6.only is a seurat object
B6.only <- readRDS('/Volumes/Samsung_T5/NKT_Analysis/RDS Files/B6NKT.rds')
Idents(B6.only) <- 'strain'

B6NKTinfo <- getEnrichedGO(B6.only, C5_df)

B6NKTinfo$NegLogPAdj <- -log10(B6NKTinfo$p.adjust)
B6NKTinfo$Cluster <- B6NKTinfo$Cluster - 1
B6NKTinfo$Cluster <- paste0('Cluster_', B6NKTinfo$Cluster)
pval.info <- as_tibble(B6NKTinfo) %>%
  group_by(Cluster) %>%
  top_n(n = -5, wt = p.adjust) %>%
  dplyr::select(-c(Description, GeneRatio, BgRatio, pvalue, qvalue, geneID, Count, p.adjust)) %>%
  drop_na() %>%
  ungroup() %>%
  pivot_wider(names_from = Cluster, values_from = NegLogPAdj, values_fill = list(NegLogPAdj = 1)) %>%
  #add_column(Cluster_6 = 1, .after = 'Cluster_5') %>% this will add columns if necessary
  transform(ID = str_replace(ID, "^GO_", "")) %>%
  column_to_rownames(var="ID") %>%
  as.matrix()
colnames(pval.info) <- gsub("_", " ", colnames(pval.info))
rownames(pval.info) <- str_replace_all(rownames(pval.info), '_', ' ')

col_fun = circlize::colorRamp2(c(2, 20), c("white", "red")) #determine range based on enrichment p-values in the list

lgd = Legend(col_fun = col_fun, title = 'padj',
             at = c(2, 5, 10, 15, 20), labels = expression(10^-2, 10^-5, 10^-10, 10^-15, 10^-20))
draw(lgd)

pval.heatmap <- Heatmap(pval.info, col = col_fun, cluster_columns = FALSE, show_row_dend = FALSE,
                        column_order = order(as.numeric(gsub("Cluster ", "", colnames(pval.info)))), 
                        cluster_rows = FALSE, row_names_gp = gpar(fontsize=13.5), border = TRUE,
                        width = unit(8, "cm"), height = unit(18, "cm"),
                        column_names_gp = gpar(fontsize=13.5))
draw(pval.heatmap, show_heatmap_legend = FALSE)