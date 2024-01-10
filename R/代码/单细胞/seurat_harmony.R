library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
# install_github('immunogenomics/harmony')
library(harmony)
library(ggplot2)
library(devtools)

dir<-paste0('00.data/',dir('00.data'))
scRNAlist<-list()
for (i in 1:length(dir)){
  ## 1. 读取数据
  # # barcodes/genes/matrix
  # Read10X()
  # h5 如果报错缺少hdf5r，运行install.packages("hdf5r")
  counts<-Read10X_h5(filename = dir[i])
  ## 2. 构建seurat对象
  pname<-strsplit(basename(dir),'.h5')[[i]]
  scRNAlist[[i]]<-CreateSeuratObject(counts = counts, project = pname,min.cells = 3, min.features = 200)
  names(scRNAlist)[i]<-pname
  scRNAlist[[i]][['tissue_type']]<-strsplit(strsplit(basename(dir),'.h5')[[i]],'_')[[1]][2]
  scRNAlist[[i]][['percent.mt']]<-PercentageFeatureSet(scRNAlist[[i]],pattern = '^MT-')
  scRNAlist[[i]][['percent.ribo']]<-PercentageFeatureSet(scRNAlist[[i]],pattern = '^RP[SL]')
}

## 3. QC
dir.create("01.QC")

theme_set1<-theme(axis.title.x = element_blank())
plot_features<-c('nFeature_RNA','nCount_RNA','percent.mt','percent.ribo')
for(i in names(scRNAlist)){
  bfqc_vln<-VlnPlot(scRNAlist[[i]],group.by = 'tissue_type',features = plot_features,ncol = 4,pt.size = 0,raster=FALSE)
  ggsave(paste0('01.QC/bfqc_vln_',i,'.pdf'),bfqc_vln,width = 10,height = 7)
  
  p1 <- FeatureScatter(scRNAlist[[i]], group.by = 'tissue_type', feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE) + NoLegend()
  p2 <- FeatureScatter(scRNAlist[[i]], group.by = 'tissue_type', feature1 = "nFeature_RNA", feature2 = "percent.mt",raster=FALSE) + NoLegend()
  p3 <- FeatureScatter(scRNAlist[[i]], group.by = 'tissue_type', feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
  bfqc_fs<-wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)
  ggsave(paste0('01.QC/bfqc_fs_',i,'.pdf'),bfqc_fs,width = 8,height = 8)
}

names(scRNAlist)
scRNAlist[["GSM6508438_cirrhotic"]]<-subset(scRNAlist[["GSM6508438_cirrhotic"]], subset = nFeature_RNA < 6000 & nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 20)
scRNAlist[["GSM6508439_HCC"]]<-subset(scRNAlist[["GSM6508439_HCC"]], subset = nFeature_RNA < 6000 & nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 15)
scRNAlist[["GSM6508440_cirrhotic"]]<-subset(scRNAlist[["GSM6508440_cirrhotic"]], subset = nFeature_RNA < 6000 & nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 15)
scRNAlist[["GSM6508441_HCC"]]<-subset(scRNAlist[["GSM6508441_HCC"]], subset = nFeature_RNA < 5000 & nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 15)


##wenzhang codes
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  print(unique(x@meta.data$orig.ident))
  # x <- subset(x, subset = percent.mt < 40 )#& nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA < 40000 )#&
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE, selection.method = "vst", nfeatures = 3000)
})



# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)

# scale and compute pca for RPCA integration method
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# perform RPCA integration of all samples
data.anchors <- FindIntegrationAnchors(object.list = scRNAlist, reduction = "rpca", anchor.features = features) # 
rm(scRNAlist)
rm(counts)
rm(bfqc_fs)
rm(bfqc_vln)
rm(p1)
rm(p2)
rm(p3)
# create an 'integrated' data assay
scRNA <- IntegrateData(anchorset = data.anchors) #, dims = 1:50



DefaultAssay(scRNA) <- "integrated"

plots = list()
for(i in seq_along(plot_features)){
  plots[[i]] = VlnPlot(scRNA,pt.size = 0,group.by = 'tissue_type', features = plot_features[i],raster=FALSE) + theme_set1 + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow=2)
ggsave(paste0('01.QC/afqc_vln_',i,'.pdf'),violin,width = 8,height = 8)


# Run the standard workflow for visualization and clustering
dir.create('02.Cluster')
scRNA <- ScaleData(scRNA, vars.to.regress = NULL,verbose = FALSE)##vars.to.regress = c(n_counts,pct)
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA, ndims = 50)
ndim <- 15
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:ndim)
scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:ndim)
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims = 1:ndim,k.param = 50)

scRNA <- FindClusters(scRNA,resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6))#, resolution = 0.5
saveRDS(scRNA,file = '01.QC/scRNA_rpca.rds')
# scRNA <- readRDS('01.QC/scRNA_rpca.rds')

DimPlot(scRNA, label = T,raster=FALSE)
p1<-DimPlot(scRNA,group.by = "group",label = F,shuffle = T,raster=FALSE)
ggsave('02.Cluster/cluster_cirr_hcc.pdf',p1,width = 8,height = 8)

for(i in colnames(scRNA@meta.data)){
  if(is.factor(scRNA@meta.data[,i])){
    levels(scRNA@meta.data[,i])<-levels(as.factor(as.numeric(levels(scRNA@meta.data[,i]))))
    p2<-DimPlot(scRNA, group.by = i,label = T,raster=FALSE)
    ggsave(paste0('02.Cluster/',i,'.pdf'),p2,width = 8,height = 8)
  } else{
    if(is.numeric(scRNA@meta.data[,i])){
      p2 <- FeaturePlot(object = scRNA, features = i,raster=FALSE)
      ggsave(paste0('02.Cluster/',i,'.pdf'),p2,width = 8,height = 8)
    } else{
      p2 <- DimPlot(scRNA,group.by = i,label = F,shuffle = T,raster=FALSE)
      ggsave(paste0('02.Cluster/',i,'.pdf'),p2,width = 8,height = 8)
    }
  }
}


##### 4 . annotation
dir.create('03.Annotation')
markerGenes <- c(#'CD74','ZEB2',##CD 4+
  'THEMIS','PRKCH','CBLB',##CD 8+
  'IGHG1','IGKC',##Plasma cell
  'IGFBP7','PDGFRA',#Mesenchymal cell
  'CFH','CYP2E1','TF','APOE',##hepatic cell
  'COL1A1','DCN',#Fibroblast
  'CD163','CD68',##Macrophage
  'TWIST2','EPCAM','BACE2',#Epithelial cell
  'FLT1','LDB2','VWF',##Endothelial cell
  'ANXA4','CFTR')##Ductal cell	Cholangiocyte
basemarker <- c('KRT18','KRT8','KRT14','KRT4','KRT19',
                'EPCAM','VWF','PECAM1',
                'COL11A1','ADGRE5',
                'CD36','RBP1','ACTA2',
                'CX3CR1','FGFBP2',
                'CD14','CD163','CD68',
                'SDC1','MZB1',
                'IGHG1','MS4A1','CD79A','CD19',
                'TPSB2','TPSAB1',
                'CD8A','CD3E','CD3D','PTPRC')

DefaultAssay(scRNA) <- "RNA"
Idents(scRNA)<-scRNA$integrated_snn_res.0.5
FeaturePlot(scRNA,features = markerGenes,raster=FALSE)
FeaturePlot(object = scRNA, features = c("HGF","COL1A1","COLEC11","MYH11"), max.cutoff = "q99", pt.size = 0.5,order = T,raster=FALSE)
DotPlot(scRNA, features = markerGenes) + RotatedAxis()
DotPlot(scRNA, features = markerGenes, cluster.idents = T,cols = c('white','red'),group.by = 'Cell.types') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  


# get the marker genes for each cluster
Idents(scRNA) <- scRNA$copykat.pred
markers <- FindAllMarkers(scRNA,only.pos = T, logfc.threshold = 0.25,max.cells.per.ident = 1000,)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing = T),]
markers$pct.1.2 <- markers$pct.1-markers$pct.2
write.csv(markers,'./03.Annotation/findallmarkers.csv')

markers <- FindMarkers(scRNA,ident.1 = 'Malignant',ident.2 = 'non_Malignant',
                       only.pos = T, logfc.threshold = 0.25,max.cells.per.ident = 1000,)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing = T),]
markers$pct.1.2 <- markers$pct.1-markers$pct.2
write.csv(markers,'./03.Annotation/copykatmarkers.csv')



FeaturePlot(object = scRNA, features = c("HGF","COL1A1","COLEC11","MYH11"), max.cutoff = "q99", pt.size = 0.5,order = T,raster=FALSE)
FeaturePlot(object = scRNA, features = c("VCAN","DCN","RGS5","LUM"), max.cutoff = "q99", pt.size = 0.5,order = T,raster=FALSE)
FeaturePlot(object = scRNA, features = c("PDGFRA","PDGFRB","LRAT","DES"), max.cutoff = "q99", pt.size = 0.5,order = T,raster=FALSE)
# 
# 
identSel <- c('Malignant','non_Malignant')##
DimPlot(scRNA, reduction = "umap", group.by = 'tissue_type',cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
FeaturePlot(scRNA,
            features = str_split("SLC7A11	 GYS1	 NDUFS1	 NDUFA11	 NUBPL	 NCKAP1	 LRPPRC	 SLC3A2	 RPN1	 ACTN4	 ACTB	 CD2AP	 CAPZB	 DSTN	 FLNA	 FLNB	 INF2	 IQGAP1	 MYH10	 MYL6	 MYH9	 PDLIM1	 TLN1",'	 ')[[1]],
            cells = unlist(CellsByIdentities(scRNA, idents = identSel)),ncol = 5)
FeaturePlot(scRNA,features = "COL1A1",cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
FeaturePlot(scRNA,features = "DCN",cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
FeaturePlot(scRNA,features = "PRKG1",cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
FeaturePlot(scRNA,features = "PDGFRB",cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
FeaturePlot(scRNA,features = "IGFBP7",cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
# FeaturePlot(object = scRNA, features = c("percent.mt"), cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
# FeaturePlot(object = scRNA, features = c("percent.ribo"), cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
# FeaturePlot(object = scRNA, features = c("nCount_RNA"), cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
# FeaturePlot(object = scRNA, features = c("nFeature_RNA"), cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
# 
# manually assigning separately clustering cells to new clusters
# https://satijalab.org/seurat/articles/visualization_vignette.html
# plot <- DimPlot(scRNA, reduction = "umap", cells = unlist(CellsByIdentities(scRNA, idents = identSel)))
# select.cells <- CellSelector(plot = plot)
# head(select.cells)
# Idents(scRNA, cells = select.cells) <- 35
# # Idents(scRNA)
# DimPlot(scRNA,label = T, pt.size = 1)
# table((scRNA@active.ident))
# levels(scRNA@active.ident)
# scRNA@active.ident <- factor(scRNA@active.ident, levels = as.character(c(0:35)))
# DimPlot(scRNA,label = T, pt.size = 1)
# scRNA[["Cell.clusters.new"]] <- Idents(object = scRNA)
# 
# # get cluster markers
# scRNA <- SetIdent(scRNA, value = "Cell.clusters.new")
# DimPlot(scRNA,label = T, pt.size = 1)
# DefaultAssay(scRNA) <- "RNA"
# DotPlot(scRNA, features = markerGenes) + RotatedAxis()
# table(scRNA@active.ident)

scRNA <- RenameIdents(object = scRNA,
                      '0' = 'Hepatocytes',
                      '1' = 'Hepatocytes',
                      '9' = 'Hepatocytes',
                      '10' = 'Hepatocytes',
                      '11' = 'Hepatocytes',
                      '12' = 'Hepatocytes',
                      '13' = 'Hepatocytes',
                      '15' = 'Hepatocytes',
                      
                      '2' = 'Endothelial',
                      
                      '3' = 'T Cell',
                      '16' = 'T Cell',
                      '8' = 'T Cell',
                      
                      '4' = 'Cholangiocytes',
                      
                      '5' = 'Mesenchymal',
                      
                      '6' = 'Fibroblast',
                      
                      '7' = 'Plasma cell',
                      
                      '14' = 'Macrophage')
scRNA[['Cell.types']] <- Idents(object = scRNA)

DimPlot(scRNA, reduction = "umap",group.by = "Cell.types", label = T,raster=FALSE)
Idents(scRNA) <- 'integrated_snn_res.0.5'
DotPlot(scRNA, features = markerGenes) + RotatedAxis()

##### .5  difference




###########harmony#################


# dir<-paste0('00.data/',dir('00.data'))
# scRNAlist<-list()
# for (i in 1:length(dir)){
#   ## 1. 读取数据
#   # # barcodes/genes/matrix
#   # Read10X()
#   # h5 如果报错缺少hdf5r，运行install.packages("hdf5r")
#   counts<-Read10X_h5(filename = dir[i])
#   ## 2. 构建seurat对象
#   pname<-strsplit(basename(dir),'.h5')[[i]]
#   scRNAlist[[i]]<-CreateSeuratObject(counts = counts, project = pname,min.cells = 3, min.features = 200)
#   scRNAlist[[i]][['tissue_type']]<-strsplit(strsplit(basename(dir),'.h5')[[i]],'_')[[1]][2]
#   scRNAlist[[i]][['percent.mt']]<-PercentageFeatureSet(scRNAlist[[i]],pattern = '^MT-')
#   scRNAlist[[i]][['percent.ribo']]<-PercentageFeatureSet(scRNAlist[[i]],pattern = '^RP[SL]')
# }
# ## 2. 整合样本
# scRNA<-merge(x = scRNAlist[[1]],y = scRNAlist[-1])
# 
# ## 3. QC
# dir.create("01.QC")
# 
# theme_set1<-theme(axis.title.x = element_blank())
# plot_features<-c('nFeature_RNA','nCount_RNA','percent.mt','percent.ribo')
# 
# bfqc_vln<-VlnPlot(scRNA,group.by = 'tissue_type',features = plot_features,ncol = 4,raster=FALSE)
# ggsave('01.QC/bfqc_vln.pdf',bfqc_vln,width = 10,height = 7)
# 
# p1 <- FeatureScatter(scRNA, group.by = 'tissue_type', feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE) + NoLegend()
# p2 <- FeatureScatter(scRNA, group.by = 'tissue_type', feature1 = "nFeature_RNA", feature2 = "percent.mt",raster=FALSE) + NoLegend()
# p3 <- FeatureScatter(scRNA, group.by = 'tissue_type', feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
# bfqc_fs<-wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)
# ggsave('01.QC/bfqc_fs.pdf',bfqc_fs,width = 8,height = 8)
# 
# scRNA<-subset(scRNA, subset = nFeature_RNA < 7500 & nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 15)
# 
# plots = list()
# for(i in seq_along(plot_features)){
#   plots[[i]] = VlnPlot(scRNA,pt.size = 0,group.by = 'tissue_type', features = plot_features[i],raster=FALSE) + theme_set1 + NoLegend()}
# violin <- wrap_plots(plots = plots, nrow=2)
# ggsave('01.QC/afqc_vln.pdf',violin,width = 8,height = 8)
# saveRDS(scRNA,file = '00.data/scRNAsubset.rds')
# 
# 
# scRNA <- NormalizeData(scRNA,normalization.method = "LogNormalize", scale.factor = 10000)
# scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
# scRNA <- ScaleData(scRNA, features = VariableFeatures(object = scRNA), vars.to.regress = c("nCount_RNA", "percent.mt"))
# scRNA <- RunPCA(scRNA, verbose = F)
# 
# # DimPlot(scRNA,reduction = 'pca',group.by = 'tissue_type')
# scRNA <- RunHarmony(scRNA,group.by.vars = 'orig.ident', max.iter.harmony = 20)##, max.iter.harmony = 20
# ElbowPlot(scRNA, ndims = 50)
# 
# scRNA <- FindNeighbors(scRNA, dims = 1:30,reduction = 'harmony',k.param = 20)##k.param = 20
# num <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6)
# for (i in num) {
#   scRNA <- FindClusters(scRNA, resolution = i,reduction = 'harmony',)
# }
# scRNA <- RunUMAP(scRNA,dims = 1:30,reduction = 'harmony',n.neighbors = 30)##,n.neighbors = 30
# DimPlot(scRNA,reduction = 'umap',group.by = 'tissue_type',split.by = 'orig.ident',raster=FALSE)##orig.ident
# DimPlot(scRNA,reduction = 'umap',group.by = 'tissue_type',split.by = 'orig.ident',ncol = 2,raster=FALSE)##orig.ident
# DimPlot(scRNA,reduction = 'umap',group.by = 'RNA_snn_res.0.3',label=T,raster=FALSE)
# 



