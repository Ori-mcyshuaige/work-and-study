library(openxlsx)
library(stringr)
library(ggplot2)
library(glmnet)
library(caret)



## 100个随机数种子
set.seed(1234)
seeds<-c(56,sample(c(1:55,57:999),99))



raw<-read.xlsx('data_all_rf.xlsx',rowNames = T)
data<-raw#[, -c(ncol(raw),ncol(raw)-1)]
data$Group<- unlist(sapply(rownames(data),function(i){strsplit(i,'_')[[1]][1]}))
data$Group<- unlist(sapply(data$Group,function(i){ifelse(i=='C','Control',i)}))
data$Group<-factor(data$Group,levels = unique(data$Group))


if (!dir.exists('Lasso-logistic')){dir.create('Lasso-logistic')}
splitdata<-data.frame(samples=rownames(data),row.names = rownames(data))
coefdata<-data.frame(samples=colnames(data)[-ncol(data)],row.names = colnames(data)[-ncol(data)])


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
  traindata<- data[sid,]
  testdata<- data[-c(sid),]
  
  #用训练集的均值和标准差来对训练集，测试集，以及后面的新样本进行标准化
  preProc<-preProcess(traindata[,-ncol(traindata)],method = c('center','scale'))
  traindata<-predict(preProc,traindata)
  testdata<-predict(preProc,testdata)
  alldata<-predict(preProc,news)[,-ncol(news)]
  
  splitdata[,paste0('seed_',seed)]<-'test'
  splitdata[sid,paste0('seed_',seed)]<-'train'
  
  
  
  ### lasso模型
  set.seed(56)
  cv.fit<-cv.glmnet(as.matrix(traindata[,-ncol(traindata)]),traindata$Group,family="binomial",type.measure='class',nfolds=10,standardize=F)
  ##变量筛选
  gt.coef <- coef(cv.fit$glmnet.fit, s = cv.fit$lambda.1se)
  
  coefdata[,paste0('seed_',seed)]<-data.frame(coefficience=gt.coef[-1],row.names =rownames(gt.coef)[-1])[rownames(coefdata),1]
  
  
  ### 预测结果
  predictdata[,paste0('seed_',seed)]<-as.numeric(unlist(predict(cv.fit,newx = as.matrix(alldata),type='response')))
}

write.xlsx(splitdata,paste0('Lasso-logistic/split.xlsx'))
write.xlsx(coefdata,paste0('Lasso-logistic/Coef.xlsx'))
write.xlsx(predictdata,paste0('Lasso-logistic/predict.xlsx'))



