library(openxlsx)
library(stringr)
library(ggplot2)
library(caret)
library(Cairo)
library(ggpubr)
library(pals)###jet颜色配色
library(ROCR)  
library(pROC)
library(ggrepel)###解决标签重叠问题







raw<-read.xlsx("data.xlsx",sheet = 1,startRow = 1,colNames = TRUE,rowNames = TRUE,check.names = T)
#colnames(raw)<-unlist(sapply(colnames(raw), function(i){paste0('`',i,'`')}))
#raw$`2',3'-Cyclic CMP `
# colnames(raw)<-unlist(sapply(colnames(raw), function(i){if(length(grep('\\d',substring(i,1,1)))==1){paste0('X',i)}
#   else if(length(grep('\\d',substring(i,1,1)))==0){i}}))
data<-raw
# data<-data.frame(log(data))
# data<-data.frame(t(scale(t(data))))
# data<-data.frame(scale(data))
data$Group<- unlist(sapply(rownames(data),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))
data$Group<-factor(data$Group,levels = unique(data$Group))
set.seed(1234)
smin<-min(unlist(sapply(unique(data$Group),function(i){sum(data$Group==i)})))
sid<-c(unlist(sapply(unique(data$Group),function(i){sample(rownames(data)[data$Group==i],3/4*smin)})))
sid<-unlist(sapply(sid, function(i){which(rownames(data)==i)}))
# sid<-sample(nrownames(data),3/4*nrow(data))
traindata<- data[sid,]
testdata<- data[-c(sid),]

# #用训练集的均值和标准差来对训练集，测试集，以及后面的新样本进行标准化
# preProc<-preProcess(traindata[,-ncol(traindata)],method = c('center','scale'))
# traindata<-predict(preProc,traindata)
# testdata<-predict(preProc,testdata)


set.seed(1234)
if (!dir.exists('SVM-REF')){dir.create('SVM-REF')}
###SVM-REF模型
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
# 执行SVM-RFE算法
set.seed(1234)
results <- rfe(traindata[, -ncol(traindata)], 
               traindata[, ncol(traindata)], 
               sizes = c(1:ncol(traindata)), 
               rfeControl = control,
               method = "svmRadial")

if (!dir.exists(paste0('SVM-REF/','0 Model RData'))){dir.create(paste0('SVM-REF/','0 Model RData'))}
save(results,file=paste0('SVM-REF/','0 Model RData/','SVM.RData'))
write.csv(traindata,paste0('SVM-REF/','0 Model RData/','traindata.csv'),row.names =T)
write.csv(testdata,paste0('SVM-REF/','0 Model RData/','testdata.csv'),row.names=T)


bios<-results$optVariables

if (!dir.exists(paste0('SVM-REF/','1 Model'))){dir.create(paste0('SVM-REF/','1 Model'))}
write.xlsx( results[["results"]],file=paste0('SVM-REF/','1 Model','/results.xlsx'),row.names = F)
write.csv( results$variables,file=paste0('SVM-REF/','1 Model','/all.csv'),row.names = F)


plotdata<-results$results
plotdata$name<-''
plotdata$name[1:length(bios)]<-unlist(sapply(bios,function(i){paste0(str_sub(i,1,20),'......')}))
p <- ggplot(data = plotdata, mapping = aes(x = Variables,y = Accuracy))+
  geom_point(aes(color=name),size=16)+
  geom_line(size=3)+
  geom_text_repel(
    data = plotdata[plotdata$name !='',], 
    mapping = aes(x = Variables,y = Accuracy,label = name,color=name), 
    size = 16,
    point.padding = unit(1, 'lines'),
    min.segment.length = unit(0.1, "lines"), 
    segment.size = 1.2,
    # segment.color = 'green', 
    force=7,box.padding=0.2,
    show.legend = F
  ) +
  theme_bw()+
  theme(
    legend.text = element_text(face = "bold", color = "black",size = ifelse(length(bios)<=10,28,56/length(bios))),
                        
    legend.title = element_text(face = "bold", color = "black",size = 36),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'black', size = nrow(traindata)*0.03),
    axis.ticks = element_line(size = 1.2,linetype=2),
    axis.ticks.length = unit(0.3,'cm'),
    axis.text.y= element_text(face = "bold", color = data$col, size = 36),
    axis.text.x= element_text(face = "bold", color = 'black', size = 36),
    axis.title = element_text(face = "bold", color = "black", size = 56)
  )
ggsave(file=paste0('SVM-REF/','1 Model/','Accuracy.pdf'),p,width = ifelse(length(bios)<=10,30,3*length(bios)),
       height = ifelse(length(bios)<=10,30,3*length(bios)),limitsize = FALSE)
ggsave(file=paste0('SVM-REF/','1 Model/','Accuracy.png'),p,width = min(50,ifelse(length(bios)<=10,30,3*length(bios))),
       height = min(50,ifelse(length(bios)<=10,30,3*length(bios))),dpi=300,limitsize = FALSE)


biodata<-traindata[,c(bios,'Group')]
if (!dir.exists(paste0('SVM-REF/','2 biomarkers'))){dir.create(paste0('SVM-REF/','2 biomarkers'))}
write.xlsx(biodata,file=paste0('SVM-REF/','2 biomarkers','/biomarkers.xlsx'),row.names = T)


# # 结果分析
# print(results)
# # 列出选择的变量集
# predictors(results)
# # 绘制结果
# plot(results, type=c("g", "o"))


###验证集ROC曲线
if (length(bios)>1) {
  set.seed(1234)
  ###用biomarkers重新绘制新模型
  set.seed(1234)
  newresults <- rfe(traindata[,bios], 
                 traindata[, ncol(traindata)], 
                 sizes = c(1:length(bios)), 
                 rfeControl = control,
                 method = "svmRadial")
  save(newresults,file=paste0('SVM-REF/','0 Model RData/','Biomarkers.RData'))
  
  #ROC曲线绘制
  if (!dir.exists(paste0('SVM-REF/','3 ROC curve'))){dir.create(paste0('SVM-REF/','3 ROC curve'))}
  getclassification<-function(data,output){
    classification<-data.frame(y=predict(newresults,newdata = as.matrix(data[,bios]),type='class'))
    colnames(classification)[1]<-'prediction'
    classification$true<-data$Group
    write.csv(classification,file=output)
  }
  # #训练集
  traindatarocs<-roc(as.numeric(traindata$Group)-1,as.numeric(predict(newresults,newdata = as.matrix(traindata[,bios]),type='response'))-1,ci=T,percent=TRUE)
  getclassification(traindata,output=paste0('SVM-REF/','3 ROC curve/','traindata.csv'))
  # #验证集
  testdatarocs<-roc(as.numeric(testdata$Group)-1,as.numeric(predict(newresults,newdata = as.matrix(testdata[,bios]),type='response'))-1,ci=T,percent=TRUE)
  getclassification(testdata,output=paste0('SVM-REF/','3 ROC curve/','testdata.csv'))
  
  
  png(paste0('SVM-REF/','3 ROC curve/','ROC plot.png'))
  plot.roc(traindatarocs,print.thres=F,print.auc=F,col='red',reuse.auc=F)
  plot.roc(testdatarocs,print.thres=F,print.auc=F,col='blue',reuse.auc=F,add = T)
  legendtest<-c('DATA: AUC(95%CI) %',
                paste0('traindata: ',sprintf("%.2f (%.2f-%.2f) %%", traindatarocs$auc, traindatarocs$ci[1],traindatarocs$ci[3])),
                paste0('testdata: ',sprintf("%.2f (%.2f-%.2f) %%", testdatarocs$auc, testdatarocs$ci[1],testdatarocs$ci[3]))
  )
  legendcol<-c('black','red','blue')
  legend('bottomright',
         legendtest,col=legendcol,text.col = legendcol,bty = 'n')
  dev.off()
 
}
file.copy("SVM-ref筛选biomarker.R",
          paste0('SVM-REF/','0 Model RData/','SVM-ref筛选biomarker.R'),
          overwrite = T)





