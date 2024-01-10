rm(list = ls())
options(stringsAsFactors = F)
gc()
library(ggpubr)
library(DESeq2)
library(ggplot2)
library(limma)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(CBNplot)
library(parallel)
#data1 =data.table::fread('./GSE203005_RAW/GSM6149650_Steady_state_pancreas_Flt3_YFPpositive_Macrophages_1_gene_counts.txt.txt.gz',data.table = F)
# files = list.files('./GSE203005_RAW/')
# dir ='./GSE203005_RAW/'
# 
# file_list = lapply(files,function(x){
#   file =paste0(dir,x)
#   data =data.table::fread(file,data.table = F)
#   data =data[,c(1,7)]
#   col_name <- strsplit(x,"_")[[1]][1]
#   colnames(data)[2] <-  col_name
#   return(data)
# })
# alldata <- do.call(cbind,file_list)
# duplicates <- duplicated(names(alldata))
# alldata <- alldata[, !duplicates]
# rownames(alldata) <- alldata[,1]
# alldata <- alldata[,-1]
# data1 <- alldata[,c(1:16)]

raw <- read.table('00.data/__MACOSX/result/gene_counts/._gene_counts.xls')
raw <- read.csv('01.节律与差异节律/全部基因logCPM.csv',row.names = 1)
data <- read.xlsx('06.TF/新建 Microsoft Excel 工作表.xlsx',sheet = 1,rowNames = T,colNames = T,check.names = F)
data <- read.csv('01.节律与差异节律/difrhygene_log2cpm+1.csv',row.names = 1)
data1 <- raw[rownames(data),]
logFoldChange=1
adjustP=0.05
for(i in colnames(data1)){
  data1[,i] <- as.numeric(data1[,i])
}
#构建表型文件
coldata <- data.frame(Sample = factor(c(rep("Ctrl",23), rep("DM",24)), 
                                      levels = c("Ctrl", "DM")),
                      row.names = colnames(data1))
dds <- DESeqDataSetFromMatrix(countData = data1 , 
                              colData = coldata, 
                              design = ~ Sample ) 

dds1 <- DESeq(dds)
v = varianceStabilizingTransformation(dds1, blind=FALSE)
vsted = assay(v)
res = results(dds1,  contrast = c('Sample', 'DM', 'Ctrl'),)
sig <- subset(res, pvalue<0.05)
incSample = rownames(subset(coldata, Sample=="DM"))


vsted <- data1[,grep('DM',colnames(data1))]
ensembl <- clusterProfiler::bitr(rownames(vsted), fromType="SYMBOL", toType="ENSEMBL", OrgDb=org.Mm.eg.db)
ensembl <- ensembl[match(unique(ensembl$SYMBOL),ensembl$SYMBOL),]
vsted <- vsted[ensembl$SYMBOL,]
rownames(vsted) <- ensembl$ENSEMBL

cand.entrez = clusterProfiler::bitr(rownames(vsted), fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)$ENTREZID
pwaykegg <- clusterProfiler::enrichKEGG(cand.entrez,organism="mmu", pvalueCutoff=1, qvalueCutoff=1)
pwaykegg@keytype <- "ENTREZID"
pwaykegg <- setReadable(pwaykegg, OrgDb=org.Mm.eg.db)
pwaykegg = enrichplot::pairwise_termsim(pwaykegg)
kegg <- data.frame(pwaykegg)
write.csv(kegg,file = 'kegg.csv')
pdf.options(reset = T)
pdf('cbnplot.pdf',width = 7,height = 7)
bngeneplot(results = pwaykegg,
                exp = vsted1,
                # expSample = incSample,
                pathNum =16, R = 50, showDir = F,
                convertSymbol = T, 
                expRow = "ENSEMBL",
                orgDb = org.Mm.eg.db,
                strThresh = 0.5,nStrength = 50)
dev.off()

