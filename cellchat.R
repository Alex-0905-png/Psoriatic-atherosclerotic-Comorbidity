library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(cowplot)
library(tidyverse)
library(CellChat)

af<-scRNA1
af$type <-Idents(af)
identity <- subset(af@meta.data, select = "celltype")
cellchat = createCellChat(object = af,meta =identity, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB = CellChatDB.human


showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use = CellChatDB
cellchat@netP$pathways
cellchat@DB = CellChatDB.use

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(object=cellchat,raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


par(mfrow = c(1,2), xpd=TRUE) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, color.use = mycolor,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, color.use = mycolor,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()



netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(2,3), remove.isolate = FALSE)











