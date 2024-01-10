#! /usr/bin/env Rscript
rm(list = ls())
options(stringsAsFactors = F) 
gc()
library(future)
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(clustree)
library(SingleR)
library(MAST)
library(SingleR)
library(Seurat)
library(tidyverse)
library(ggrepel)

#harmony整合结果评估
load('/home/xuqingpeng/240/00.data/single_data/scRNA_harmony.Rdata')
scRNA<- FindNeighbors(scRNA1, reduction = "harmony",dims = 1:30)
rm(scRNA1)
gc()
scRNA <- FindClusters(scRNA, resolution = 0.8, verbose = FALSE)

#####计算差异
# 多线程运行减少时间
plan()
plan("multiprocess", workers = 35)
plan()
diff.mast <- FindAllMarkers(scRNA,test.use = 'MAST',only.pos = TRUE)
save(diff.mast,file = '/home/xuqingpeng/240/05.cell type/diff.mast.RData')
##singleR 分类
refdata <- celldex::HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, 
                      stringsAsFactors = F)
scRNA@meta.data$singleR = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'singleR'] <- celltype$celltype[i]}
dir.create('05.cell type')
p0 = DimPlot(scRNA, label=T, label.size=5, reduction='umap')
p1 = DimPlot(scRNA, group.by="singleR", label=T, label.size=5, reduction='umap')
p3 <- p0|p1
ggsave("/home/xuqingpeng/240/singR.pdf", p3, width = 20, height = 10,dpi = 1000)


###人工注释
celltype_marker=c(
  "EPCAM",#上皮细胞 epithelial
  "PLVAP",#内皮细胞 endothelial
  "LUM","DCN",#成纤维细胞 fibroblasts
  "KIT",##mast cell
  'PLD4',##dendritic cell 
  "CD163",#macrophage
  "MS4A1",#B细胞
  "TNFRSF17",#浆细胞 plasma cell
  'KLRD1',#NK细胞
  "CD8A","GZMH","GZMM","NKG7"#T细胞
  )
p4 <- VlnPlot(scRNA,features = celltype_marker,pt.size = 0,ncol = 2)
ggsave("/home/xuqingpeng/240/marker_vinplot.pdf",p4,width = 44,height = 45,units = "cm")

(n=length(unique(scRNA@meta.data$seurat_clusters)))
cell_type=data.frame(ClusterID=0:(n-1),
                    celltype='unkown')  #构建数据框

## 判断亚群ID属于那类细胞
cell_type[celltype$ClusterID %in% c(2,8,14,17,22,30,26),2]='Epithelial'
cell_type[celltype$ClusterID %in% c(24,5),2]='Endothelial'
cell_type[celltype$ClusterID %in% c(7,12,15,27,28,29,33),2]='Fibroblasts' 
cell_type[celltype$ClusterID %in% c(11,34),2]='mast cell' 
cell_type[celltype$ClusterID %in% c(19),2]='dendritic cell'
cell_type[celltype$ClusterID %in% c(4,10),2]='macrophag'
cell_type[celltype$ClusterID %in% c(3,16,18,20,25,31,32),2]='plasma cell'
cell_type[celltype$ClusterID %in% c(9),2]='B cell'
cell_type[celltype$ClusterID %in% c(6),2]='NK cell'
cell_type[celltype$ClusterID %in% c(0,1,13,21,23),2]='T cell'

## 重新赋值
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(cell_type)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == cell_type$ClusterID[i]),'celltype'] <- cell_type$celltype[i]}
##绘图
mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","celltype"))
cell_type_med <- mydata %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
source('/home/xuqingpeng/240/umap_plot.R',print.eval = T)
p5 <- umap_plot(mydata,celltype)
ggsave("/home/xuqingpeng/240/UMAP_celltype.pdf", p5, width = 12, height = 8,dpi = 1000)

#p6 <- VlnPlot(scRNA,features = celltype_marker,pt.size = 0,ncol = 2,group.by="celltype")
#ggsave("/home/xuqingpeng/240/vlnplot_marker.pdf", p6, width = 44, height = 45,dpi = 1000,units = "cm")

mycolor <- c("lightgrey", "blue","seagreen2")#设置颜色
p7 <- FeaturePlot(scRNA, features = celltype_marker, 
                  cols = mycolor,reduction = "umap",pt.size = 1.5,ncol = 3)
ggsave("/home/xuqingpeng/240/marker_xx.pdf", p7, width = 44, height = 45,dpi = 1000,units = "cm")
