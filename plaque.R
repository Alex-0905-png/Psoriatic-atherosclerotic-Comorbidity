library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

dir = c('2As_vessel/','3As_vessel/','4As_vessel/')
names(dir) = c('2As_vessel','3As_vessel','4As_vessel')

counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB_genes, rownames(scRNA@assays$RNA))
HB_genes <- rownames(scRNA@assays$RNA)[HB_m]
HB_genes <- HB_genes[!is.na(HB_genes)]
scRNA[["HB_percent"]] <- PercentageFeatureSet(scRNA, features=HB_genes)

As <- subset(scRNA,
             subset = nFeature_RNA > 300 & nFeature_RNA < 7000 &
               mt_percent < 10 &
               HB_percent < 3 &
               nCount_RNA < quantile(nCount_RNA,0.97) & nCount_RNA > 1000)


dir = c('1Aspso_vessel/','2Aspso_vessel/')
names(dir) = c('1Aspso_vessel','2Aspso_vessel')

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



As$batch <- 'As'
Aspso$batch <- 'Aspso'

#harmony
library(harmony)
combined <- merge(As, y = Aspso, add.cell.ids = c("As", "Aspso"))
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

cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10,
              11,12,13,14,15,16,17,18,19,20,
              21,22,23,24,25,26,27,28,29,30,
              31,32,33),
  celltype = c("T_cells","T_cells","T_cells","T_cells","Macrophages","Fibroblasts",
               "Macrophages","Macrophages","T_cells","T_cells","Macrophages",
               "Macrophages","Endothelial_cells","Macrophages","T_cells","Macrophages",
               "Myeloid_cells","Macrophages","Fibroblasts","Macrophages","Macrophages",
               "Macrophages","Fibroblasts","T_cells","B_cells","Myeloid_cells",
               "B_cells","Macrophages","Endothelial_cells","Macrophages","Macrophages",
               "B_cells","Macrophages","Macrophages")
)

cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10,
              11,12,13,14,15,16,17,18,19,20,
              21,22,23,24,25,26,27,28,29,30,
              31,32,33),
  celltype = c("T_cells","T_cells","T_cells","T_cells","IL17RA-Macrophages","Fibroblasts",
               "IL17RA+Macrophages","IL17RA-Macrophages","T_cells","T_cells","IL17RA-Macrophages",
               "IL17RA+Macrophages","Endothelial_cells","IL17RA+Macrophages","T_cells","IL17RA+Macrophages",
               "Myeloid_cells","IL17RA+Macrophages","Fibroblasts","IL17RA+Macrophages","IL17RA+Macrophages",
               "IL17RA+Macrophages","Fibroblasts","T_cells","B_cells","Myeloid_cells",
               "B_cells","IL17RA+Macrophages","Endothelial_cells","IL17RA+Macrophages","IL17RA+Macrophages",
               "B_cells","IL17RA+Macrophages","IL17RA+Macrophages")
)

scRNA1@meta.data$celltype ="NA"
for(i in 1:nrow(cluster_types)) {
  cluster_id <- cluster_types$cluster[i]
  cell_type <- cluster_types$celltype[i]
  scRNA1$celltype[scRNA1$seurat_clusters == cluster_id] <- cell_type
}


DimPlot(scRNA1, reduction = "umap",group.by = "celltype",
        #cols = color1,
        pt.size = 0.5,
        label = F,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

makers<-c("CD79A","IGHM","MS4A1",
          "VWF","PECAM1","EDN1",
          "COL1A1","COL1A2","COL3A1",
          "LGMN","CD14","MSR1",
          "CTSG","TPSAB1","MS4A2",
          "CD2","IL7R","IL32")


library(MySeuratWrappers)

outFile="Rplots.pdf"
pdf(file=outFile,width=9,height=6)
VlnPlot(scRNA1, features = makers,pt.size=0,group.by = "celltype",cols = color1,
        stacked=T,direction = "horizontal")
dev.off()








library(reshape2)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)


df <- table(scRNA1@meta.data$celltype,scRNA1@meta.data$batch) %>% melt()
colnames(df) <- c("Cluster","Sample","Number")
df$Cluster <- factor(df$Cluster)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_color <- c("#FFBE7A","#FA7F6F",
                  "#63b2ee","#76da91","#7CAE00",
                  "#f89588","#9192ab","#C77CFF","#FFA500")
ggplot(data = df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )










