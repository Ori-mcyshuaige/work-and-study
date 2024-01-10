library(randomForest)
library(stringr)
library(ggpubr)
library(pals)###jet颜色配色
#library(doParallel) # 并行处理的包
library(Cairo)
library(openxlsx)
library(pROC)
#library(MUVR) # MUVR包




## 100个随机数种子
set.seed(8)
seeds1<-sample(c(1:9999),300)
set.seed(24)
seeds2<-sample(c(1:9999),1)


raw<-read.xlsx('Differentially expressed metabolites-Merge-2.xlsx',rowNames = T,sheet = 5,check.names = T,sep.names = " ")
data<-raw
for(ii in 1:ncol(data)){
  result<-sapply(as.numeric(data[,ii]),function(i){(i-mean(as.numeric(data[,ii])))/sd(as.numeric(data[,ii]))})
  data[,ii]<-result
}
data$Group<- unlist(sapply(rownames(data),function(i){strsplit(i,'_')[[1]][1]}))
data$Group<- unlist(sapply(data$Group,function(i){ifelse(i=='C','Control',i)}))
data$Group<-factor(data$Group,levels = unique(data$Group))


if (!dir.exists('300Random-forest')){dir.create('300Random-forest')}
splitdata<-data.frame(samples=rownames(data),row.names = rownames(data))
accuracydata<-data.frame(samples=colnames(data)[-ncol(data)],row.names = colnames(data)[-ncol(data)])
ginidata<-data.frame(samples=colnames(data)[-ncol(data)],row.names = colnames(data)[-ncol(data)])


# ###对离群样本的处理
# news<-read.xlsx('data_all_yz.xlsx',rowNames = T)
# news$Group<- unlist(sapply(rownames(news),function(i){strsplit(i,'_')[[1]][1]}))
# news$Group<- unlist(sapply(news$Group,function(i){ifelse(i=='C','Control',i)}))
# news<-rbind(data,news)
predictdata<-data.frame(samples=rownames(data),row.names = rownames(data))


for (seed1 in seeds1) {
  ### 划分训练集，测试集
  set.seed(seed1)
  smin<-unlist(sapply(unique(data$Group),function(i){sum(data$Group==i)}))
  names(smin)<-unique(data$Group)
  sid<-c(sample(rownames(data)[data$Group==names(smin)[1]],4/5*smin[1]),sample(rownames(data)[data$Group==names(smin)[2]],4/5*smin[2]))
  sid<-unlist(sapply(sid, function(i){which(rownames(data)==i)}))
  traindata<- data[sid,]
  testdata<- data[-c(sid),]
  
  splitdata[,paste0('seed_',seed1)]<-'test'
  splitdata[sid,paste0('seed_',seed1)]<-'train'
  
  ### 随机森林模型
  for(seed2 in seeds2){
  set.seed(seed2)
  data.rf = randomForest(Group ~ ., data=traindata, importance=TRUE, proximity=TRUE,ntree=500)
  importdata<-as.data.frame(data.rf$importance)
  if (!dir.exists('300Random-forest/seed')){dir.create('300Random-forest/seed')}
  write.xlsx(importdata,paste0('300Random-forest/seed/seed_',seed2,"_",seed1,'.xlsx'),rowNames = T,overwrite = T)
  
  accuracydata[,paste0('seed_',seed2,"_",seed1)]<-importdata[rownames(accuracydata),3]
  ginidata[,paste0('seed_',seed2,"_",seed1)]<-importdata[rownames(accuracydata),4]
  
  ### 预测结果
  predictdata[,paste0('seed_',seed2,"_",seed1)]<-predict(data.rf,newdata = data[,-ncol(data)],type='prob')[,2]
  }
}

accuracydata<-accuracydata[,-1]
ginidata<-ginidata[,-1]

gini<-NULL
accuracy<-NULL
for(i in colnames(ginidata)){
  ginidata<-ginidata[order(ginidata[,i],decreasing = T),]
  if(all(is.na(match(c("Linoleamide","Oleamide","Adrenoylethanolamide"),rownames(ginidata)[1:6])))){
    gini<-c(gini,i)
  }
  accuracydata<-accuracydata[order(accuracydata[,i],decreasing = T),]
  if(all(is.na(match(c("Linoleamide","Oleamide","Adrenoylethanolamide"),rownames(accuracydata)[1:6])))){
    accuracy<-c(accuracy,i)
  }
}
result<-data.frame(gini=gini,accuracy=c(accuracy,NA,NA,NA,NA))
write.xlsx(result,"result.xlsx")

write.xlsx(splitdata,paste0('300Random-forest/split.xlsx'))
write.xlsx(accuracydata,paste0('300Random-forest/Accuracy.xlsx'),rowNames=T,overwrite = T)
write.xlsx(ginidata,paste0('300Random-forest/Gini.xlsx'),rowNames=T,overwrite = T)
write.xlsx(predictdata,paste0('300Random-forest/predict.xlsx'))


