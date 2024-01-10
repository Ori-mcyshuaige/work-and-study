library(openxlsx)


###行为样本，列为物质
data1<-read.xlsx("data.xlsx",sheet = 1,colNames = T,rowNames = T,check.names = F,sep.names = " ")
data2<-read.xlsx("data.xlsx",sheet = 2,colNames = T,rowNames = T,check.names = F,sep.names = " ")
data2<-data2[rownames(data1),]
zdata1<-as.data.frame(matrix(ncol = ncol(data1),nrow = nrow(data1)))
colnames(zdata1)<-colnames(data1)
rownames(zdata1)<-rownames(data1)
zdata2<-as.data.frame(matrix(ncol = ncol(data2),nrow = nrow(data2)))
colnames(zdata2)<-colnames(data2)
rownames(zdata2)<-rownames(data2)
##对列做z-score
for(ii in colnames(data1)){
  result<-sapply(data1[,ii],function(i){(i-mean(data1[,ii]))/sd(data1[,ii])})
  zdata1[,ii]<-result
}
for(ii in colnames(data2)){
  result<-sapply(data2[,ii],function(i){(i-mean(data2[,ii]))/sd(data2[,ii])})
  zdata2[,ii]<-result
}
##对行做z-score
# for(ii in rownames(data)){
#       result<-sapply(as.numeric(data[ii,]),function(i){(i-mean(as.numeric(data[ii,])))/sd(as.numeric(data[ii,]))})
#       data[ii,]<-result
#     }

zdata<-cbind(zdata1,zdata2)
write.xlsx(zdata,"zdata.xlsx",rowNames = T)
