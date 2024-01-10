library(tidyverse)
library(tidyr)
library(dplyr)
library(reprex)
library(stringr)

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(robustbase)


oridata <- read.csv('./00.data/GSE98793.csv',header = T)
oridata <- oridata[,c(1,grep('no',oridata[54677,]))]
coldata <- data.frame('group'=unlist(oridata[54676,2:129]))
rownames(coldata) <- colnames(oridata)[2:129]
rawdata <- oridata[1:54675,]
# rawdata <- read.table('miRNA_exp_deseq2_input.txt',header = T)
# coldata <- read.table('coldata.txt',header = T,row.names = 'sample')
##过滤表达量低的基因
# counts_value<-apply(rawdata[,-1],1,table)
# number0<-c()
# for(i in counts_value){
#   number0<-c(number0,max(i))-
# }
# filterdata<-rawdata[number0<3,]
filterdata<-rawdata
##去重,取最大值
inputdata<-data.frame(matrix(ncol = (ncol(filterdata)-1), nrow = length(unique(filterdata[,1]))))
colnames(inputdata)<-colnames(filterdata)[-1]
rownames(inputdata)<-unique(filterdata[,1])

for(i in unique(filterdata[,1])){
  idata<-filterdata[filterdata[,1]==i,-1]
  inputdata[i,]<-colMaxs(as.matrix(idata))
}
# inputdata<-as.matrix(filterdata[,-1])
# rownames(inputdata)<-filterdata[,1]

inputdata<-as.matrix(inputdata)
coldata[,1]<-factor(coldata[,1])
##构建dds矩阵
dds <- DESeqDataSetFromMatrix(countData = inputdata, colData = coldata, design = ~ group)
dds$group <- relevel(dds$group, ref = 'Control')
##差异分析
dds <- DESeq(dds)
resname<-paste0('group_',levels(dds$group)[2],'_vs_',levels(dds$group)[1])
res <- results(dds, alpha = 0.05, name = resname)
res_df <- as.data.frame(res)
# res_df$gene_id = rownames(res_df)
statdf<-cbind(inputdata[rownames(res_df),],res_df)
write.csv(statdf,file = paste0('stat_',resname,'.csv'))
write.csv(res_df,file = paste0(resname,'.csv'))
dif_mirna<-res_df[!is.na(res_df$padj),]
dif_mirna<-dif_mirna[dif_mirna$padj<0.05&abs(dif_mirna$log2FoldChange)>1,]
dif_mirna<-cbind(inputdata[rownames(dif_mirna),],dif_mirna)
write.csv(dif_mirna,paste0('dif_',resname,'.csv'))
