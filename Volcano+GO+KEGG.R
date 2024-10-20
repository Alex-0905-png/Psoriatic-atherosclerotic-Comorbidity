library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)

dge.cluster <- FindMarkers(scRNA1,ident.1 = c(1,4),ident.2 = c(3,19))

#GO
colnames(deg)
logFC=0.5
P.Value = 0.05
type1 = (deg$p_val< P.Value)&(deg$avg_log2FC< -logFC)
type2 = (deg$p_val< P.Value)&(deg$avg_log2FC> -logFC)
deg$Group2 = ifelse(type1,"Down",ifelse(type2,"Up","Not-Sig"))
table(deg$Group2)

diff.genes<-rownames(subset(deg,Group2!='Not-Sig'))


diff.df <- bitr(diff.genes,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

go.diff2 <- enrichGO(gene = diff.df$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     pAdjustMethod = 'BH',
                     pvalueCutoff =0.01,
                     qvalueCutoff = 0.05,
                     ont="all",
                     readable =T)

deg$SYMBOL<-diff.genes
Enrichment<-deg%>%
  mutate(SYMBOL=rownames(.))%>%
  inner_join(diff.df, by = "SYMBOL")%>%
  dplyr::select(SYMBOL,ENTREZID,avg_log2FC)

go <- data.frame(Category = go.diff2$ONTOLOGY,
                 ID = go.diff2$ID,
                 Term = go.diff2$Description,
                 Genes = gsub("/", ", ", go.diff2$geneID),
                 adj_pval = go.diff2$p.adjust)

genelist <- data.frame(ID = Enrichment$SYMBOL,
                       logFC = Enrichment$avg_log2FC)

circ <- circle_dat(go, genelist)

GOCircle(circ,nsub=10)

#Volcano
res<-dge.cluster
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'avg_log2FC',
                y = 'p_val',
                drawConnectors = F)

#KEGG
res<-diff.genes
library(org.Hs.eg.db)

entrezid_all = mapIds(x = org.Hs.eg.db, 
                      keys = res, 
                      keytype = "SYMBOL", 
                      column = "ENTREZID") 

entrezid_all = na.omit(entrezid_all)  
entrezid_all = data.frame(entrezid_all) 
head(entrezid_all)
library("clusterProfiler")
library("org.Hs.eg.db")

KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], 
                         keyType = "kegg",
                         organism= "human",  
                         pAdjustMethod = "fdr", 
                         pvalueCutoff = 1,  
                         qvalueCutoff =1)  
KEGG_enrich  = data.frame(KEGG_enrich)

display_number = 12
KEGG_enrich = as.data.frame(KEGG_enrich)[1:display_number[1], ]

kk = as.data.frame(KEGG_enrich)
rownames(kk) = 1:nrow(kk)
kk$order=factor(rev(as.integer(rownames(kk))),labels = rev(kk$Description))

ggplot(kk,aes(y=order,x=Count,fill=pvalue))+  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "KEGG Pathways Enrichment",  #设置标题、x轴和Y轴名称
       x = "Gene number",
       y = "Pathway")


ggplot(kk,aes(y=order,x=GeneRatio))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(pvalue,size="Count"),
       x="Gene Ratio",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
