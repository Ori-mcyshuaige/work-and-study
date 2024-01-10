library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
###准备想要计算的gene list

scRNA<-readRDS('scRNA.rds')

genelist<-readxl::read_xlsx('PAN_genelist.xlsx')

genelist<-as.list(genelist)

scRNA <- AddModuleScore(object = scRNA, features = genelist,ctrl = 100,name = 'PAN_Score')
scRNA[['PAN_stage']] <- ifelse(scRNA@meta.data[,'PAN_Score']>mean(scRNA@meta.data[,'PAN_Score']),'High','Low')
Idents(scRNA) <- 'PAN_Score'
DimPlot(scRNA,reduction = 'umap',label = TRUE,pt.size = .2)


