library(openxlsx)
library(stringr)
library(ggplot2)
library(glmnet)
library(Cairo)
library(ggpubr)
library(pals)###jet颜色配色
library(ROCR)  
library(pROC)
library(caret)



scatterplot<-function(i,data,df.c){
  x <<- as.character(df.c$x[i])
  y <<- as.character(df.c$y[i])
  library(ggpubr)
  if (which(colnames(data)==x)<which(colnames(data)==y)) {
    p <- ggscatter(data, x = x,y = y,shape = 20,size = 5,
                   color = "Group",
                   palette = jet(length(levels(data$Group))+2)[2:(length(levels(data$Group))+1)][as.numeric(unique(data$Group))],add.params = list(color = "black"),
                   add = "reg.line",cor.method='spearman',
                   conf.int = T, cor.coef = T,show.legend.text = F)+
      geom_rug(aes(colour=Group))+
      #scale_color_manual(values=jet(length(levels(data$Group))+2)[2:(length(levels(data$Group))+1)][as.numeric(unique(data$Group))])+
      theme_bw() +
      theme(legend.position=c(.925,.5),
            #title =  element_text(face = "bold",color = "black",size = 24),
            axis.title = element_text(face = "bold",color = "black",size = 16),
            axis.ticks= element_blank(),
            #axis.text = element_blank(),
            #axis.text.x = element_text(face = "bold",color = "black",size = 16),
            legend.key = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = 'black', size = 1))
  }else if (which(colnames(data)==x)>which(colnames(data)==y)) {
    datax<-data[,c(x,"Group")]
    datax$meta<-x
    colnames(datax)[1]<-'exep'
    datay<-data[,c(y,"Group")]
    datay$meta<-y
    colnames(datay)[1]<-'exep'
    data<-rbind(datax,datay)
    #colnames(data)[1]<-'exep'
    ###均值
    descdata1<-as.data.frame(aggregate(data[,'exep'],list(meta=data[,'meta'],Group=data[,'Group']),mean))
    ##标准误
    descdata2<-as.data.frame(aggregate(data[,'exep'],list(meta=data[,'meta'],Group=data[,'Group']),function(i){sd(i)/sqrt(length(i))}))
    ##95%置信空间
    descdata3<-as.data.frame(aggregate(data[,'exep'],list(meta=data[,'meta'],Group=data[,'Group']),function(i){(sd(i)/sqrt(length(i)))*qt(.95/2 + .5, length(i)-1)}))
    descdata<-cbind.data.frame(descdata1,descdata2[,3],descdata3[,3])
    colnames(descdata)[3:5]<-c('mean','se','ci')
    data<-merge(data,descdata,by=c('meta','Group'))
    
    p<-ggplot(descdata,aes(x=Group,y=mean,fill=Group,colour=meta,group=meta,linetype=meta))+
      geom_bar(stat = 'identity',position=position_dodge(0.4),width=.3)+
      #geom_point(data=data,aes(y=exep),position=position_dodge(0.4),size=1.2,alpha=0.5)+
      geom_line(position=position_dodge(0.4),size=.56)+
      geom_errorbar(aes(ymin=mean-ci,ymax=mean+ci),position=position_dodge(0.4),size=.56,width=.36)+
      annotate('text',label=unique(descdata$meta)[1],x=2,y=max(descdata$mean)*.95,size=12*.95,face='blod',color=jet(2)[1])+
      annotate('text',label=unique(descdata$meta)[2],x=2,y=max(descdata$mean)*.85,size=12*.95,face='blod',color=jet(2)[2])+
      scale_fill_manual(values = jet(length(levels(data$Group))+2)[2:(length(levels(data$Group))+1)][as.numeric(unique(data$Group))])+
      scale_color_manual(values = jet(2))+
      scale_shape_manual(values = c(21,22))+
      theme_bw()+
      theme(legend.position='none',
            axis.title = element_blank(),
            axis.ticks= element_blank(),
            #axis.text.y = element_blank(),
            axis.text.x = element_text(face = "bold",color = "black",size = 16),
            legend.key = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = 'black', size = 1))
  }else{
    data$meta<-data[,y]
    descdata<-as.data.frame(aggregate(data[,'meta'],list(Group=data[,'Group']),mean))
    
    p<-ggplot(data,aes(x=as.numeric(Group)-.3,y=meta,fill=Group,color=Group,group=Group))+
      geom_violin()+
      geom_boxplot(width=.1,notch=T,fill='white')+
      geom_dotplot(aes(x=as.numeric(Group)+.3,group=Group),binaxis='y',stackdir = 'center',method='histodot')+
      geom_rug(sides='left')+
      annotate('text',label=y,x=2,y=max(data[,y])*.95,size=12*.95,face='blod')+
      scale_fill_manual(values=jet(length(levels(data$Group))+2)[2:(length(levels(data$Group))+1)][as.numeric(unique(data$Group))])+
      scale_color_manual(values=jet(length(levels(data$Group))+2)[2:(length(levels(data$Group))+1)][as.numeric(unique(data$Group))])+
      scale_x_continuous(breaks=1:nlevels(data$Group),labels =levels(data$Group))+
      theme_bw()+
      theme(legend.position=c(.925,.5),
            #title =  element_text(face = "bold",color = "black",size = 24),
            axis.title = element_blank(),
            axis.ticks= element_blank(),
            #axis.text.y = element_blank(),
            axis.text.x = element_text(face = "bold",color = "black",size = 16),
            legend.key = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = 'black', size = 1))
  }
  
  
  return(p)
}
Features<-function(data,output=''){
  library(grid)
  library(gridExtra)
  df.c <- data.frame(
    x = rep(colnames(data)[-ncol(data)], each = ncol(data)-1), 
    y = rep(colnames(data)[-ncol(data)], ncol(data)-1)
  )
  l.pics <- lapply(
    1:nrow(df.c), scatterplot, 
    df.c = df.c,data = data
  )
  names(l.pics) <- paste0('p', df.c$x, df.c$y)
  p.grid <- arrangeGrob(grobs = l.pics, nrow = ncol(data)-1)
  ggsave(paste0(output,'Features plot.pdf'),p.grid,width=56,height=56,limitsize = FALSE)
  # ggsave('Features plot.svg',p,width=56,height=56,limitsize = FALSE)
  # CairoPNG(file="Features plot.png",width=(ncol(data)+1)*56,height=(ncol(data)+1)*56)
  # grid.draw(p.grid)
  # dev.off()
  return(p.grid)
}


importanceplot<-function(data,output=''){
  
  #options(stringsAsFactors = F)
  require(scales)
  getcolorname_scale_color_gradient2<-function(x,low='blue',high='red',mid='#999999',space="Lab",midpoint=0){
    require(scales)
    get_palette<-div_gradient_pal(low,mid, high, space)
    
    get_mid_rescaler<-function(x=x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
      rescale_mid(x, to, from, mid = midpoint)
    }
    return(get_palette(get_mid_rescaler(x)))
  }
  options(stringsAsFactors = T)
  data$col<-getcolorname_scale_color_gradient2(data$coefficience)

  p<-ggplot(data,aes(x=biomarkers,y=coefficience,color=coefficience,fill=coefficience))+
    geom_bar(aes(x=biomarkers,y=coefficience,color=coefficience,fill=coefficience),stat = "identity",width = 0.9,alpha = 0.9) +
    scale_color_gradient2(low='blue',high='red',mid='#999999',name = 'coefficience')+
    scale_fill_gradient2(low='blue',high='red',mid='#999999',name = 'coefficience')+
    #title('Regression Coefficience')+
    ylab('coefficience')+
    coord_flip()+
    theme_bw() +
    theme(
      legend.text = element_text(face = "bold", color = "black", size = 28,family = 'Times'),
      legend.title = element_text(face = "bold", color = "black", size = 36,family = 'Times'),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'black', size = 1),
      axis.ticks = element_line(size = 1.2,linetype=2),
      axis.ticks.length = unit(0.3,'cm'),
      axis.text.y= element_text(face = "bold", color = data$col, size = 28),
      axis.text.x= element_text(face = "bold", color = 'black', size = 28),
      axis.title= element_text(face = "bold", color = "black", size = 36),
      title = element_text(face = "bold", color = "black", size = 36)
    )
  ggsave(
    paste0(output, 'varImpPlot.jpg'), p, 
    width = 20, height = max(0.2*nrow(data),8), units = 'in', dpi = 300
  )
  ggsave(
    paste0(output, 'varImpPlot.pdf'), p, 
    width = 20, height = max(0.2*nrow(data),8), units = 'in', dpi = 300
  )
}


# rocplot<-function(data,output){
#   #tp<-sum(data$prediction==data$true)
#   pred <- prediction(data$prediction,data$true)   #预测值(0.5二分类之前的预测值)和真实值     
#   auc<-performance(pred,'auc')@y.values        #AUC值  
#   perf <- performance(pred,'tpr','fpr')  
#   #plot(perf) 
#   auc_ci=roc(data$true,data$prediction,ci=TRUE)
#   data<-data.frame(x=perf@x.values[[1]],y=perf@y.values[[1]])
#   p<-ggplot(data,aes(x,y))+
#     geom_line(color='red')+
#     theme_bw()+
#     theme(
#       panel.grid.major = element_blank(), 
#       panel.grid.minor = element_blank(),
#       
#     )
# }












raw<-read.xlsx('RF.xlsx',sheet = 1,rowNames = T)
#colnames(raw)<-unlist(sapply(colnames(raw), function(i){paste0('`',i,'`')}))
#raw$`2',3'-Cyclic CMP `
# colnames(raw)<-unlist(sapply(colnames(raw), function(i){if(length(grep('\\d',substring(i,1,1)))==1){paste0('X',i)}
#   else if(length(grep('\\d',substring(i,1,1)))==0){i}}))
data<-raw#[, -c(ncol(raw),ncol(raw)-1)]
for(ii in 1:ncol(data)){
  result<-sapply(as.numeric(data[,ii]),function(i){(i-mean(as.numeric(data[,ii])))/sd(as.numeric(data[,ii]))})
  data[,ii]<-result
}
# data<-data.frame(log(data))
# data<-data.frame(t(scale(t(data))))
# data<-data.frame(scale(data))
data$Group<- unlist(sapply(rownames(data),function(i){strsplit(i,'_')[[1]][1]}))
data$Group<- unlist(sapply(data$Group,function(i){ifelse(i=='C','Control',i)}))
data$Group<-factor(data$Group,levels = unique(data$Group))
set.seed(1234)
smin<-min(unlist(sapply(unique(data$Group),function(i){sum(data$Group==i)})))
sid<-c(unlist(sapply(unique(data$Group),function(i){sample(rownames(data)[data$Group==i],3/4*smin)})))
sid<-unlist(sapply(sid, function(i){which(rownames(data)==i)}))
# sid<-sample(nrownames(data),3/4*nrow(data))
traindata<- data[sid,]
testdata<- data[-c(sid),]

#用训练集的均值和标准差来对训练集，测试集，以及后面的新样本进行标准化
preProc<-preProcess(traindata[,-ncol(traindata)],method = c('center','scale'))
traindata<-predict(preProc,traindata)
testdata<-predict(preProc,testdata)

###LASSO
if (!dir.exists('Lasso-logistic')){dir.create('Lasso-logistic')}
set.seed(1234)
cv.fit<-cv.glmnet(as.matrix(traindata[,-ncol(traindata)]),traindata$Group,family="binomial",type.measure='class',nfolds=10,standardize=F)
if (!dir.exists(paste0('Lasso-logistic/','1 Model'))){dir.create(paste0('Lasso-logistic/','1 Model'))}
png(paste0('Lasso-logistic/','1 Model','/independent variable.png'))
plot(cv.fit)
dev.off()
pdf(paste0('Lasso-logistic/','1 Model','/independent variable.pdf'))
plot(cv.fit)
dev.off()
if (!dir.exists(paste0('Lasso-logistic/','0 Model RData'))){dir.create(paste0('Lasso-logistic/','0 Model RData'))}
save(cv.fit,file=paste0('Lasso-logistic/','0 Model RData/','Lasso-logistic.RData'))
write.csv(traindata,paste0('Lasso-logistic/','0 Model RData/','traindata.csv'),row.names =T)
write.csv(testdata,paste0('Lasso-logistic/','0 Model RData/','testdata.csv'),row.names=T)

##变量筛选
gt.coef <- coef(cv.fit$glmnet.fit, s = cv.fit$lambda.1se)
#gt.coef[which(gt.coef != 0)] 
bios<-rownames(gt.coef)[which(gt.coef != 0)][-1]
bioplot<-data.frame(coefficience=gt.coef[which(gt.coef != 0)][-1],biomarkers=bios,row.names =bios )
# bioplot$coefficiences<-unlist(sapply(bioplot$coefficience,function(i){if(i>0){log(i)}else{-log(abs(i))}}))
bioplot<-bioplot[order(bioplot$coefficience),]
bioplot$biomarkers<-factor(bioplot$biomarkers,levels = bioplot$biomarkers)
importanceplot(bioplot,output = paste0('Lasso-logistic/','1 Model/'))
write.xlsx(bioplot,file=paste0('Lasso-logistic/','1 Model','/coefficience.csv'),row.names = T)

biodata<-traindata[,c(bios,'Group')]
if (!dir.exists(paste0('Lasso-logistic/','2 biomarkers'))){dir.create(paste0('Lasso-logistic/','2 biomarkers'))}
write.xlsx(biodata,file=paste0('Lasso-logistic/','2 biomarkers','/biomarkers.xlsx'),row.names = T)
# p<-Features(biodata,output = paste0('Lasso-logistic/','2 biomarkers/'))
# CairoPNG(file="Features plot.png",width=(ncol(traindata)+1)*56,height=(ncol(traindata)+1)*56)
# grid.draw(p)
# dev.off()
# file.copy("Features plot.png",paste0('Lasso-logistic/','2 biomarkers/','Features plot.png'),overwrite=T)
# file.remove("Features plot.png")



###模型作图 判断是否符合logistic
if (!dir.exists(paste0('Lasso-logistic/','3 logistic regression'))){dir.create(paste0('Lasso-logistic/','3 logistic regression'))}
predata<-data.frame(y=predict(cv.fit,newx = as.matrix(traindata[,-ncol(traindata)]),type='response'),Group=traindata$Group)
colnames(predata)[1]<-'prediction'

predata<-predata[order(predata$prediction),]

predata$sample<-seq(0,1,length.out = nrow(traindata))
p <- ggplot(data = predata, mapping = aes(x = sample,y = prediction,color=Group))+
  geom_point(size=9)+
  theme_bw()+
  theme(
    legend.text = element_text(face = "bold", color = "black", size = nrow(traindata)*0.36),
    legend.title = element_text(face = "bold", color = "black", size = nrow(traindata)*0.56),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'black', size = nrow(traindata)*0.03),
    axis.ticks = element_blank(),
    axis.text= element_blank(),
    axis.title = element_text(face = "bold", color = "black", size = nrow(traindata)*0.56)
  )
ggsave(file=paste0('Lasso-logistic/','3 logistic regression/','logistic regression.pdf'),p,width = 30,height = 30)
ggsave(file=paste0('Lasso-logistic/','3 logistic regression/','logistic regression.png'),p,width = 30,height = 30,dpi=300)



rocdata<-data.frame(predict(cv.fit,newx = as.matrix(traindata[,-ncol(traindata)]),type='class'))
colnames(rocdata)<-'prediction'
rocdata$true<-traindata$Group
write.csv(rocdata,file=paste0('Lasso-logistic/','3 logistic regression/','logistic regression.csv'))



###验证集ROC曲线
if (length(bios)>1) {
  set.seed(1234)
  ###用biomarkers重新绘制新模型
  cv.fittest<-cv.glmnet(as.matrix(traindata[,bios]),traindata$Group,family="binomial",type.measure='class',nfolds=10,standardize=F)
  save(cv.fittest,file=paste0('Lasso-logistic/','0 Model RData/','Biomarkers.RData'))
  gt.coeftest <- coef(cv.fittest$glmnet.fit, s = cv.fittest$lambda.1se)
  
  #ROC曲线绘制
  if (!dir.exists(paste0('Lasso-logistic/','4 ROC curve'))){dir.create(paste0('Lasso-logistic/','4 ROC curve'))}
  getclassification<-function(data,output){
    classification<-data.frame(y=predict(cv.fittest,newx = as.matrix(data[,bios]),type='class'))
    colnames(classification)[1]<-'prediction'
    classification$true<-data$Group
    classification$prob<-unlist(predict(cv.fittest,newx = as.matrix(data[,bios]),type='response'))
    write.csv(classification,file=output)
  }
  # #训练集
  traindatarocs<-roc(as.numeric(traindata$Group)-1,predict(cv.fittest,newx = as.matrix(traindata[,bios]),type='response'),ci=T,percent=TRUE)
  getclassification(traindata,output=paste0('Lasso-logistic/','4 ROC curve/','traindata.csv'))
  # #验证集
  testdatarocs<-roc(as.numeric(testdata$Group)-1,predict(cv.fittest,newx = as.matrix(testdata[,bios]),type='response'),ci=T,percent=TRUE)
  getclassification(testdata,output=paste0('Lasso-logistic/','4 ROC curve/','testdata.csv'))
  # #得到新样本，观察模型预测新样本能力
  # 女性离群样本
  # liqundataraw<-read.xlsx('data_yz.xlsx',rowNames = T)
  # liqundataraw[,colnames(traindata)[-ncol(traindata)]]<-predict(preProc,liqundataraw[,colnames(traindata)[-ncol(traindata)]])#用训练集的均值和标准差来对新样本进行标准化
  # liqundata<-liqundataraw[, bios]
  # liqundata$Group<- unlist(sapply(rownames(liqundata),function(i){strsplit(i,'_')[[1]][1]}))
  # liqundata$Group<- unlist(sapply(liqundata$Group,function(i){ifelse(i=='C','Control',i)}))
  # liqundata$Group<-factor(liqundata$Group,levels = unique(liqundata$Group))
  # 
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   liqundata$Group<-factor(liqundata$Group,levels = unique(liqundata$Group))
  #   liqundataroc<-roc(as.numeric(liqundata$Group)-1,predict(cv.fittest,newx = as.matrix(liqundata[,bios]),type='response'),ci=T,percent=TRUE)
  }
  # getclassification(liqundata,output=paste0('Lasso-logistic/','4 ROC curve/','离群.csv'))
  
  
  png(paste0('Lasso-logistic/','4 ROC curve/','ROC plot.png'))
  plot.roc(traindatarocs,print.thres=F,print.auc=F,col='red',reuse.auc=F)
  plot.roc(testdatarocs,print.thres=F,print.auc=F,col='blue',reuse.auc=F,add = T)
  legendtest<-c('DATA: AUC(95%CI) %',
                paste0('traindata: ',sprintf("%.2f (%.2f-%.2f) %%", traindatarocs$auc, traindatarocs$ci[1],traindatarocs$ci[3])),
                paste0('testdata: ',sprintf("%.2f (%.2f-%.2f) %%", testdatarocs$auc, testdatarocs$ci[1],testdatarocs$ci[3]))
  )
  legendcol<-c('black','red','blue')
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
  #   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
  #   legendcol<-c(legendcol,'green')
  # }
  legend('bottomright',
         legendtest,col=legendcol,text.col = legendcol,bty = 'n')
  dev.off()

file.copy("Lasso-logistic回归筛选biomarker.R",
          paste0('Lasso-logistic/','0 Model RData/','Lasso-logistic回归筛选biomarker.R'),
          overwrite = T)

