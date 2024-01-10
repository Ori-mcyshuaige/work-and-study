library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(rio)
library(GSEABase)


data <- read.csv('01.差异分析/GSE98793.csv',row.names = 1)
data <- data[,grep('MDD',colnames(data))]
#ssGSEA分析
geneset <- qusage::read.gmt("immune.gmt")
# genelist <- data.table::fread("./gene.csv",data.table = F)
# geneset$disulfidptosis <- genelist[[1]]
# colnames(genelist)[2] <- "type"
# genelist <- lapply(split(genelist,genelist$type), function(x){
#   dd = x$gene
#   unique(dd)
# })
data <- data[-grep('group',rownames(data)),]
for(i in colnames(data)){
  data[,i] <- as.numeric(data[,i])
}
ssgseaScore=gsva(as.matrix(data), geneset, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
write.csv(ssgseaScore,file="ssGSEA.mianyi.result.csv",row.names = T)

##boxplot
data1 <- data.frame(sample=rep(colnames(ssgseaScore),nrow(ssgseaScore)),
                    `Fraction Score`=unlist(sapply(rownames(ssgseaScore),function(x){rbind(as.data.frame(unlist(ssgseaScore[x,])))})),
                    celltype=rep(rownames(ssgseaScore),each=ncol(ssgseaScore))
                    ,check.names = F)
data1$Cluster <- 'C1'
data1$Cluster[unlist(lapply(nameC2,function(x){grep(x,data1$sample)}))] <- 'C2'
data1$celltype <- str_replace_all(data1$celltype,'\\.',' ')
#HLA
HLAgene <- c('HLA-A','HLA-B','HLA-C','HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB',
             'HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DRA','HLA-DRB4',
             'HLA-E','HLA-F','HLA-G')
data <- data[HLAgene,]
data2 <- data.frame(sample=rep(colnames(data),nrow(data)),
                    `Expression Value`=unlist(sapply(rownames(data),function(x){rbind(as.data.frame(unlist(data[x,])))})),
                    MDDGene=rep(rownames(data),each=ncol(data))
                    ,check.names = F)
data2$Cluster <- 'C1'
data2$Cluster[unlist(lapply(nameC2,function(x){grep(x,data2$sample)}))] <- 'C2'
data2$MDDGene <- str_replace_all(data2$MDDGene,'\\.',' ')

p=ggboxplot(data2, x="MDDGene", y="Expression Value", color = "Cluster",
            xlab="",
            ylab="Gene expression",
            legend.title="Cluster",
            palette = crgCluCol,
            width=0.6,
            add="point")+rotate_x_text(60)+
  stat_compare_means(aes(group=Cluster),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                     label = "p.signif")
ggsave("./MDDboxplot.pdf",p,width=11, height=7,dpi = 600)




#####################################################

data =data[,c(98:194)]
geneSets=getGmt('./step5-fig5/c2.cp.kegg.symbols.gmt', geneIdType=SymbolIdentifier())
ssgseaScore=gsva(data, geneSets, method='gsva')

#对GSVA的打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

#读取分型的结果文件
cluster=read.table('./step5-fig5/cluster.txt', header=T, sep="\t", check.names=F, row.names=1)
nameC1=row.names(cluster[cluster$Cluster=="C1",,drop=F])
nameC2=row.names(cluster[cluster$Cluster=="C2",,drop=F])
dataC1=ssgseaScore[,nameC1,drop=F]
dataC2=ssgseaScore[,nameC2,drop=F]
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
data=cbind(dataC1, dataC2)
Type=c(rep("C1",conNum), rep("C2",treatNum))