library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

##Aspso data
dir = c('Aspso/','Aspso2/')
names(dir) = c('Aspso','Aspso2/')
counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB_genes, rownames(scRNA@assays$RNA))
HB_genes <- rownames(scRNA@assays$RNA)[HB_m]
HB_genes <- HB_genes[!is.na(HB_genes)]
scRNA[["HB_percent"]] <- PercentageFeatureSet(scRNA, features=HB_genes)

Aspso <- subset(scRNA,
                subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
                  mt_percent < 10 &
                  HB_percent < 3 &
                  nCount_RNA < quantile(nCount_RNA,0.97) & nCount_RNA > 1000)

#Pso data
dir = c('Pso1/', "Pso2/")
names(dir) = c('Pso1',  'Pso2')
counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB_genes, rownames(scRNA@assays$RNA))
HB_genes <- rownames(scRNA@assays$RNA)[HB_m]
HB_genes <- HB_genes[!is.na(HB_genes)]
scRNA[["HB_percent"]] <- PercentageFeatureSet(scRNA, features=HB_genes)

Pso <- subset(scRNA,
              subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
                mt_percent < 10 &
                HB_percent < 3 &
                nCount_RNA < quantile(nCount_RNA,0.97) & nCount_RNA > 1000)

#Con data
dir = c('Con1/', "Con2/", "Con3/", "Con4/", "Con5/")
names(dir) = c('Con1',  'Con2',  'Con3',  'Con4',  'Con5')
counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB_genes, rownames(scRNA@assays$RNA))
HB_genes <- rownames(scRNA@assays$RNA)[HB_m]
HB_genes <- HB_genes[!is.na(HB_genes)]
scRNA[["HB_percent"]] <- PercentageFeatureSet(scRNA, features=HB_genes)

Con <- subset(scRNA,
              subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
                mt_percent < 10 &
                HB_percent < 3 &
                nCount_RNA < quantile(nCount_RNA,0.97) & nCount_RNA > 1000)

Aspso$batch <- 'Aspso'
Pso$batch <- 'Pso'
Con$batch <- 'Con'

#harmony
library(harmony)
combined <- merge(Aspso, y = list(Pso, Con), add.cell.ids = c("Aspso", "Pso", "Con"))
#combined<-Pso
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
head(combined@meta.data)

harmony_integrated <- RunHarmony(combined, group.by.vars = "batch", plot_convergence = TRUE)
scRNA1<-harmony_integrated
#scRNA1<-combined
pc.num=1:20
scRNA1<-FindNeighbors(scRNA1,dims = pc.num)
scRNA1<-FindClusters(scRNA1,resolution = 1.2)
scRNA1<-BuildClusterTree(scRNA1)

# UMAP
scRNA1 <- RunUMAP(scRNA1, reduction = "harmony", dims = 1:20)

#figure show
DimPlot(scRNA1, reduction = "umap",#group.by = "batch",
        #cols = color1,
        pt.size = 0.5,
        label = F,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

#maker
markers <- FindAllMarkers(object = scRNA1, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10,
              11,12,13,14,15,16,17,18,19,20,
              21,22,23,24,25,26,27,28,29,30,
              31,32,33,34,35,36,37,38,39),
  celltype = c("Keratinocytes", "Endothelial_cells", "Endothelial_cells","Fibroblasts","Mast_cells","Keratinocytes",
               "T_cells","T_cells","T_cells","Keratinocytes","Fibroblasts",
               "Fibroblasts","Fibroblasts","Keratinocytes","Fibroblasts","Macrophages",
               "T_cells","Keratinocytes","Keratinocytes","Keratinocytes","Endothelial_cells",
               "Macrophages","Endothelial_cells","Macrophages","Keratinocytes","Endothelial_cells",
               "Keratinocytes","Macrophages","T_cells","Macrophages","T_cells",
               "Fibroblasts","Mast_cells","Keratinocytes","T_cells","Macrophages",
               "Macrophages","Keratinocytes","Macrophages","Mast_cells")
)

scRNA1@meta.data$celltype ="NA"
for(i in 1:nrow(cluster_types)) {
  cluster_id <- cluster_types$cluster[i]
  cell_type <- cluster_types$celltype[i]
  scRNA1$celltype[scRNA1$seurat_clusters == cluster_id] <- cell_type
}









celltypemakers<-c("PIP","SCGB1B2P","MUCL1",
                  "CDH5","ICAM1","PECAM1",
                  "COL1A1","COL1A2","COL3A1",
                  "S100A9","KRT23","S100A7",
                  "TPSAB1","TPSB2","CTSG",
                  "LYZ","HLA-DQA1","CD68",
                  "CD2","CD3D","IL32")


DoHeatmap(scRNA1,
          features = as.character(unique(celltypemakers)),
          group.by = "celltype",
          assay = "RNA",
          group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90","#00CD00","#008B8B","#FFA500"))+
  scale_fill_gradientn(colors = c("#58539f","white","#d86967"))











