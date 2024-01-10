library(xgboost)
library(openxlsx)
library(stringr)
library(ggplot2)



## 100个随机数种子
set.seed(1234)
seeds<-c(56,sample(c(1:55,57:999),99))


raw<-read.xlsx('data_all_rf.xlsx',rowNames = T)
data<-raw#[, -c(ncol(raw),ncol(raw)-1)]
data$Group<- unlist(sapply(rownames(data),function(i){strsplit(i,'_')[[1]][1]}))
data$Group<- unlist(sapply(data$Group,function(i){ifelse(i=='C','Control',i)}))
data$Group<-factor(data$Group,levels = unique(data$Group))


if (!dir.exists('XGBoost-tree')){dir.create('XGBoost-tree')}
splitdata<-data.frame(samples=rownames(data),row.names = rownames(data))
gaindata<-data.frame(samples=colnames(data)[-ncol(data)],row.names = colnames(data)[-ncol(data)])
coverdata<-data.frame(samples=colnames(data)[-ncol(data)],row.names = colnames(data)[-ncol(data)])
frequencydata<-data.frame(samples=colnames(data)[-ncol(data)],row.names = colnames(data)[-ncol(data)])



###对离群样本的处理
news<-read.xlsx('data_all_yz.xlsx',rowNames = T)
news$Group<- unlist(sapply(rownames(news),function(i){strsplit(i,'_')[[1]][1]}))
news$Group<- unlist(sapply(news$Group,function(i){ifelse(i=='C','Control',i)}))
news<-rbind(data,news)
predictdata<-data.frame(samples=rownames(news),row.names = rownames(news))


for (seed in seeds) {
  ### 划分训练集，测试集
  set.seed(seed)
  smin<-min(unlist(sapply(unique(data$Group),function(i){sum(data$Group==i)})))
  sid<-c(unlist(sapply(unique(data$Group),function(i){sample(rownames(data)[data$Group==i],3/4*smin)})))
  sid<-unlist(sapply(sid, function(i){which(rownames(data)==i)}))
  traindataraw<- data[sid,]
  testdataraw<- data[-c(sid),]
  
  #用训练集的均值和标准差来对训练集，测试集，以及后面的新样本进行标准化
  ###转换为xgboost特定的数据模式
  traindata1<-Matrix(data.matrix(traindataraw[,c(-ncol(traindataraw))]),sparse=T)
  traindata<-xgb.DMatrix(data=traindata1,label = as.numeric(traindataraw$Group)-1)
  
  splitdata[,paste0('seed_',seed)]<-'test'
  splitdata[sid,paste0('seed_',seed)]<-'train'
  
  
  
  ### XGboost模型
  set.seed(56)
  xgboostmodel<-xgboost(data = traindata,max_depth=6, eta=0.7,  objective='binary:logistic',booster='gbtree', nround=100,seed=56)
  ##变量筛选
  bios=xgb.importance(model = xgboostmodel)
  
  gaindata[,paste0('seed_',seed)]<-sapply(rownames(gaindata), function(i){ifelse(i %in% bios$Feature,bios$Gain[bios$Feature==i],0)})
  coverdata[,paste0('seed_',seed)]<-sapply(rownames(coverdata), function(i){ifelse(i %in% bios$Feature,bios$Cover[bios$Feature==i],0)})
  frequencydata[,paste0('seed_',seed)]<-sapply(rownames(frequencydata), function(i){ifelse(i %in% bios$Feature,bios$Frequency[bios$Feature==i],0)})
  
  ### 预测结果
  alldata<-Matrix(data.matrix(news[,-ncol(news)]),sparse=T)
  predictdata[,paste0('seed_',seed)]<-predict(xgboostmodel,newdata = alldata)
  
}

write.xlsx(splitdata,paste0('XGBoost-tree/split.xlsx'))
write.xlsx(gaindata,paste0('XGBoost-tree/Gain.xlsx'))
write.xlsx(coverdata,paste0('XGBoost-tree/Cover.xlsx'))
write.xlsx(frequencydata,paste0('XGBoost-tree/Frequency.xlsx'))
write.xlsx(predictdata,paste0('XGBoost-tree/predict.xlsx'))



