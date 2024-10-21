library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(png)
library(jsonlite)


spe2 = Read10X("./matrix/")
image2 <- Read10X_Image(image.dir = file.path("./", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
#没有报错，无需转化
for (i in colnames((spe2@images$slice1@coordinates))) {
  spe2@images$slice1@coordinates[[i]] <- as.integer(spe2@images$slice1@coordinates[[i]])
}

SpatialFeaturePlot(spe2, features = "nFeature_Spatial")


SpatialFeaturePlot(spe2, features ="IL17RA")
SpatialFeaturePlot(spe2, features ="PPARG")







