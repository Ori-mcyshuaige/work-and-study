library(ggm)
library(psych)
library(openxlsx)

setwd("D:/马驰宇/桌面/何院专辑/研究院项目/戴主任项目重分析/偏相关")
data1<-read.xlsx("临床.xlsx",sheet = 4,rowNames = T)###包含控制条件
data2<-read.xlsx("蛋白.xlsx",sheet = 1,rowNames = T)
# data2<-data2[rownames(data1),]
data<-cbind(data1,data2)

wb <-createWorkbook()
cor<-matrix(nrow = length(data1)-3,ncol = length(data2))
cor<-data.frame(cor)
rownames(cor)<-colnames(data1)[-c(1:3)]
colnames(cor)<-colnames(data2)
pvalue<-cor
s<-cov(data)
for(i in colnames(data1)[-c(1:3)]){
  for(j in colnames(data2)){
    jsbl<-c(i,j)
    tjbl<-colnames(data1[,1:3])
    u<-c(jsbl,tjbl)
    r<-pcor(u,s)
    p<-pcor.test(r,length(tjbl),dim(data)[1])$pvalue
    cor[i,j]<-r
    pvalue[i,j]<-p
  }
}
p<-NULL
a<-1
qvalue<-matrix(nrow = nrow(pvalue),ncol = ncol(pvalue))
qvalue<-data.frame(qvalue)
rownames(qvalue)<-rownames(pvalue)
colnames(qvalue)<-colnames(pvalue)
for(kk in 1:length(pvalue)){
  z<-pvalue[,kk]
  p<-c(p,z)}
q<-p.adjust(p,method = "BH")
# q<-fdrtool(p,statistic="pvalue")$lfdr
for(ii in 1:ncol(pvalue)){
  for(jj in 1:nrow(pvalue)){
    qvalue[jj,ii]<-q[a]
    a<-a+1
  }
}
X1<-NULL
X2<-NULL
Corr<-NULL
Pvalue<-NULL
Qvalue<-NULL
for(m in 1:nrow(pvalue)){
  for(n in 1:ncol(pvalue)){
    if(pvalue[m,n]<0.05){
      X1<-c(X1,rownames(pvalue)[m])
      X2<-c(X2,colnames(pvalue)[n])
      Corr<-c(Corr,cor[m,n])
      Pvalue<-c(Pvalue,pvalue[m,n])
      Qvalue<-c(Qvalue,qvalue[m,n])
      select<-data.frame(X1,X2,Corr,Pvalue,Qvalue)
    }
  }
}
addWorksheet(wb,sheetName = "cor",gridLines = T)
addWorksheet(wb,sheetName = "pvalue",gridLines = T)
addWorksheet(wb,sheetName = "qvalue",gridLines = T)
addWorksheet(wb,sheetName = "select",gridLines = T)
writeDataTable(wb,sheet = 1,cor,rowNames = T)
writeDataTable(wb,sheet = 2,pvalue,rowNames = T)
writeDataTable(wb,sheet = 3,qvalue,rowNames = T)
writeDataTable(wb,sheet = 4,select,rowNames = T)
saveWorkbook(wb,overwrite = T,"DAI蛋白偏相关.xlsx")


