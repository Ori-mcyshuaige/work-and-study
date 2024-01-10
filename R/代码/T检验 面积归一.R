library(openxlsx)
setwd("C:/Users/小马很酷/Desktop/To李溢溪")
ori<-read.xlsx("NEG-Mean升降.xlsx",sheet = 2,colNames = T,rowNames = T,check.names = F,sep.names = " ")
gydata<-ori[,grep("Control|Psoriasis",colnames(ori))]

##面积归一
gydata <- sapply(gydata, function(column) {
  sum.col <- sum(column)
  column <- column / sum.col
})

gydata<-as.data.frame(gydata)
datac<-ori[,grep("A",colnames(ori))]
# datas<-gydata[,grep("Psoriasis",colnames(gydata))]
datas<-ori[,c(1:length(ori))[-c(grep("A",colnames(ori)))]]
pvalue<-NULL
meanc<-NULL
means<-NULL

##T-test  &  Mean
for(i in 1:nrow(datac)){
  var<-var.test(unlist(datac[i,],use.names = F),unlist(datas[i,],use.names = F))$p.value > 0.05
  p<-t.test(unlist(datac[i,],use.names = F), unlist(datas[i,],use.names = F), var.equal = var)$p.value
  mean1<-mean(unlist(datac[i,],use.names = F))
  mean2<-mean(unlist(datas[i,],use.names = F))
  pvalue<-c(pvalue,p)
  meanc<-c(meanc,mean1)
  means<-c(means,mean2)
}

qvalue<-p.adjust(pvalue,method = "BH")
output<-data.frame(meanNoRenal=meanc,meanRenal=means,fc=means/meanc,pvalue=pvalue,qpvalue=qvalue)
rownames(output)<-rownames(ori)
# output<-cbind(datac,datas)
# output<-cbind(ori[,c(1,2,3)],output)
# output$meanc<-meanc
# output$means<-means
# output$pvalue<-pvalue
# output$fc<-means/meanc
# output$logfc<-log2(output$fc)

write.xlsx(output,"RESULT2.xlsx",rowNames=T,overwrite = T)
