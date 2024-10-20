library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(stringr)
library(cowplot)
library(scales)
library(readr)
library(progeny)
library(gplots)
library(tibble)
library(grid)
library(rlang)

Macrophages_M1<-c( 'IL6', 'NOS2', 'TLR2',  'TLR4',
                   "CD80","CD86")


library(Hmisc)
Macrophages_M1_gene=capitalize(tolower(Macrophages_M1)) %>% list()

sce<-sub_pbmc02
sce =  AddModuleScore(object = sce,features = Macrophages_M1_gene )
colnames(sce@meta.data)
plot(sce$Cluster1,sce$Cluster2,col=sce$celltype)

All.merge=AddModuleScore(sce,
                         features = Macrophages_M1_gene,
                         name = "Inflammatory Score")
colnames(All.merge@meta.data)

VlnPlot(All.merge,features = "Inflammatory.Score1")
VlnPlot(All.merge,features = "Inflammatory.Score1",pt.size = 0,group.by = "celltype",combine = TRUE)




















