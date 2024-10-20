library(monocle)
pbmc<-scRNA1
pbmc<-sub_pbmc02
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = pbmc@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))


diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~ celltype")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.1)) ## 不要也写0.1 ，而是要写0.01。

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

HSMM <- reduceDimension(HSMM, max_components = 3,
                        num_dim = 20,
                        method = "DDRTree") 


HSMM <- orderCells(HSMM)

plot_cell_trajectory(HSMM)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

plot_cell_trajectory(HSMM, color_by = "orig.ident")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

plot_cell_trajectory(HSMM, color_by = "celltype")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())




plot_cell_trajectory(HSMM, color_by = "celltype")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

pData(HSMM)$CXCL8 <- log2(exprs(HSMM)['CXCL8',]+1)
plot_cell_trajectory(HSMM, color_by = "CXCL8") +
  scale_color_gsea()+theme_dr(xlength = 0.22, ylength = 0.22,
                              arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())


all<-HSMM

blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("PPARG")))
plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "seurat_clusters",
                  min_expr = 0.1)
plot_genes_violin(HSMM[blast_genes,],
                  grouping = "seurat_clusters",
                  min_expr = 0.1)


HSMM_expressed_genes <-  row.names(subset(fData(HSMM_myo),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("PPARG", "IL17RA")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Hours")



library(patchwork)
topN <- head(ordering_genes, 3)
topN <- c("CXCL8","IL17RA")

topN_cds <- HSMM[topN,]

p1 <-plot_genes_in_pseudotime(topN_cds, color_by = "Pseudotime")
p2 <-plot_genes_in_pseudotime(topN_cds, color_by = "celltype")

