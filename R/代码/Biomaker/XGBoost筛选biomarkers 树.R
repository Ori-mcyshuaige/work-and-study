library(xgboost)
library(openxlsx)
library(stringr)
library(ggplot2)
library(Cairo)
library(ggpubr)
library(pals)###jet颜色配色
library(ROCR)  
library(pROC)
library(Matrix)
library(DiagrammeR)





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
  getcolorname_scale_color_gradient<-function(x,low='#999999',high='red',space="Lab"){
    require(scales)
    get_palette<-seq_gradient_pal(low, high, space)
    
    get_mid_rescaler<-function(x=x, to = c(-1, 1), from = range(x, na.rm = TRUE)) {
      rescale_mid(x, to, from)
    }
    return(get_palette(get_mid_rescaler(x)))
  }
  options(stringsAsFactors = T)
  data<-data[order(data$Gain,decreasing = F),]
  data$Feature<-factor(data$Feature,levels = data$Feature)
  data$col<-getcolorname_scale_color_gradient(data$Gain,low='#FFBFBF')
  
  p<-ggplot(data,aes(x=Feature,y=Gain,color=Gain,fill=Gain))+
    geom_bar(aes(x=Feature,y=Gain,color=Gain,fill=Gain),stat = "identity",width = 0.9,alpha = 0.9) +
    scale_color_gradient(low='#FFBFBF',high='red',name = 'coefficience')+
    scale_fill_gradient(low='#FFBFBF',high='red',name = 'coefficience')+
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

cvplot<-function(data,output=''){
  train<-data[,c('iter','train_error_mean','train_error_std')]
  colnames(train)<-c('iter','mean','std')
  train$type<-'traindata'
  test<-data[,c('iter','test_error_mean','test_error_std')]
  colnames(test)<-c('iter','mean','std')
  test$type<-'testdata'
  data<-rbind(train,test)
  pd<-position_dodge(0.1)
  p<-ggplot(data,aes(x=iter,y=mean,colour=type,group=type))+
    # # geom_errorbar(aes(ymin=mean-std,ymax=mean+std),position = pd,size=0.05,width=0.05)+
    # geom_point(position = pd,size=0.3)+
    geom_line(position = pd,size=1.7)+
    ggtitle('10 fold cross validation')+
    ylab('mean(error)')+
    ylim(0,max(data$mean))+
    theme_bw() +
    theme(
      legend.text = element_text(face = "bold", color = "black", size = 28,family = 'Times'),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'black', size = 1),
      axis.ticks = element_line(size = 1.2,linetype=2),
      axis.ticks.length = unit(0.3,'cm'),
      axis.text.y= element_text(face = "bold", color = data$col, size = 28),
      axis.text.x= element_text(face = "bold", color = 'black', size = 28),
      axis.title= element_text(face = "bold", color = "black", size = 36),
      plot.title = element_text(face = "bold", color = "blue", size = 36,hjust = 0.5,vjust = 0.5)
    )
  ggsave(
    paste0(output, '10 fold cross validation.jpg'), p, 
    width = 20, height = 20, units = 'in', dpi = 300
  )
  ggsave(
    paste0(output, '10 fold cross validation.pdf'), p, 
    width = 20, height = 20, units = 'in', dpi = 300
  )
}






raw<-read.xlsx("data.xlsx",sheet = 1,startRow = 1,colNames = TRUE,rowNames = TRUE,check.names = T)
#colnames(raw)<-unlist(sapply(colnames(raw), function(i){paste0('`',i,'`')}))
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
traindataraw<- data[sid,]
testdataraw<- data[-c(sid),]
###转换为xgboost特定的数据模式
traindata1<-Matrix(data.matrix(traindataraw[,c(-ncol(traindataraw))]),sparse=T)
traindata<-xgb.DMatrix(data=traindata1,label = as.numeric(traindataraw$Group)-1)
# testdata1<-Matrix(data.matrix(testdataraw[,c(-ncol(testdataraw))]),sparse=T)
# testdata<-xgb.DMatrix(data=testdata1,label = as.numeric(testdataraw$Group)-1)


###XGBoost
if (!dir.exists('XGBoost-tree')){dir.create('XGBoost-tree')}
set.seed(1234)
xgboostmodel<-xgboost(data = traindata,max_depth=6, eta=0.7,  objective='binary:logistic',booster='gbtree', nround=100,seed=1234)
# xgboostmodelliner<-xgboost(data = traindata,max_depth=9, eta=0.56,  objective='binary:logistic',
#                            booster='gblinear', nround=100,seed=1234, callbacks = list(cb.gblinear.history()))

if (!dir.exists(paste0('XGBoost-tree/','1 Model'))){dir.create(paste0('XGBoost-tree/','1 Model'))}
gr <- xgb.plot.tree(model=xgboostmodel, trees=NULL, render=FALSE,show_node_id = T)
export_graph(gr, paste0('XGBoost-tree/','1 Model','/tree.png'), width=1200, height=10800)
export_graph(gr, paste0('XGBoost-tree/','1 Model','/tree.pdf'), width=1200, height=3600)
# # <unknown>:1919791: Invalid asm.js: Function definition doesn't match use
# # https://github.com/rich-iannone/DiagrammeRsvg/issues/7
# # 坑爹呢！export_graph调用R包DiagrammeRsvg，而DiagrammeRsvg的作者Rich 有bug，作者自己都没办法解决
# apuf<-tempfile()
# # cat(DiagrammeRsvg::export_svg(grViz(apuf)),file=apuf)
# rsvg::rsvg_pdf(apuf,"koe2.pdf",width=5500,height=5500/sqrt(2))
# rsvg::rsvg_png(gr,"koe2.png",width=5500,height=5500/sqrt(2))



if (!dir.exists(paste0('XGBoost-tree/','0 Model RData'))){dir.create(paste0('XGBoost-tree/','0 Model RData'))}
save(xgboostmodel,file=paste0('XGBoost-tree/','0 Model RData/','XGBoost-tree.RData'))
write.csv(traindataraw,paste0('XGBoost-tree/','0 Model RData/','traindata.csv'),row.names =T)
write.csv(testdataraw,paste0('XGBoost-tree/','0 Model RData/','testdata.csv'),row.names=T)


##变量筛选
bios=xgb.importance(model = xgboostmodel)
importanceplot(bios,output = paste0('XGBoost-tree/','1 Model/'))
write.xlsx(bios,file=paste0('XGBoost-tree/','1 Model','/importance.xlsx'),row.names = F)


# 10折交叉验证
if (!dir.exists(paste0('XGBoost-tree/','2 10fold cross validation'))){dir.create(paste0('XGBoost-tree/','2 10fold cross validation'))}
xgboostmodelcv<-xgb.cv(data=traindata,max_depth=6,eta=0.7,objective='binary:logistic',booster='gbtree', nround=100,seed=1234,nfold = 10, showsd=T)
itertoerror<-xgboostmodel$evaluation_log
write.csv(itertoerror,paste0('XGBoost-tree/','2 10fold cross validation/','itertoerror.csv'),row.names =F)
itertoerrorcv<-xgboostmodelcv$evaluation_log
write.csv(itertoerrorcv,paste0('XGBoost-tree/','2 10fold cross validation/','10 fold cross validation.csv'),row.names =F)
cvplot(itertoerrorcv,output = paste0('XGBoost-tree/','2 10fold cross validation/'))



###biomarker
biodata<-traindataraw[,c(bios$Feature,'Group')]
if (!dir.exists(paste0('XGBoost-tree/','3 biomarkers'))){dir.create(paste0('XGBoost-tree/','3 biomarkers'))}
write.xlsx(biodata,file=paste0('XGBoost-tree/','3 biomarkers','/biomarkers.xlsx'),row.names = T)
# p<-Features(biodata,output = paste0('XGBoost-tree/','3 biomarkers/'))
# CairoPNG(file="Features plot.png",width=(ncol(traindata)+1)*56,height=(ncol(traindata)+1)*56)
# grid.draw(p)
# dev.off()
# file.copy("Features plot.png",paste0('XGBoost-tree/','3 biomarkers/','Features plot.png'),overwrite=T)
# file.remove("Features plot.png")






###验证集ROC曲线
if (length(bios$Feature)>1) {
  set.seed(1234)
  ###转换为xgboost特定的数据模式
  traindata1<-Matrix(data.matrix(traindataraw[,bios$Feature]),sparse=T)
  newtrain<-xgb.DMatrix(data=traindata1,label = as.numeric(traindataraw$Group)-1)
  testdata1<-Matrix(data.matrix(testdataraw[,bios$Feature]),sparse=T)
  newtest<-xgb.DMatrix(data=testdata1,label = as.numeric(testdataraw$Group)-1)
  ### 重新构造模型
  xgboostmodelnew<-xgboost(data = newtrain,max_depth=6, eta=0.7,  objective='binary:logistic',booster='gbtree', nround=100,seed=1234)
  save(xgboostmodelnew,file=paste0('XGBoost-tree/','0 Model RData/','Biomarkers.RData'))
  if (!dir.exists(paste0('XGBoost-tree/','4 ROC curve'))){dir.create(paste0('XGBoost-tree/','4 ROC curve'))}
  
  
  getclassification<-function(data,raw,output){
    classification<-data.frame(y=round(predict(xgboostmodelnew,newdata = data)))
    colnames(classification)[1]<-'prediction'
    classification$true<-raw$Group
    numtogroup<-unique(traindataraw$Group)
    classification$prediction<-unlist(sapply(classification$prediction,function(i){as.character(numtogroup[[i+1]])}))
    rownames(classification)<-rownames(raw)
    classification$prob<-predict(xgboostmodelnew,newdata = data)
    write.csv(classification,file=output)
  }
  # #训练集
  traindatarocs<-roc(as.numeric(traindataraw$Group)-1,predict(xgboostmodelnew,newdata = newtrain),ci=T,percent=TRUE)
  getclassification(newtrain,traindataraw,output=paste0('XGBoost-tree/','4 ROC curve/','traindata.csv'))
  # #验证集
  testdatarocs<-roc(as.numeric(testdataraw$Group)-1,predict(xgboostmodelnew,newdata = newtest),ci=T,percent=TRUE)
  getclassification(newtest,testdataraw,output=paste0('XGBoost-tree/','4 ROC curve/','testdata.csv'))
  
  png(paste0('XGBoost-tree/','4 ROC curve/','ROC plot.png'))
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


file.copy("XGBoost筛选biomarkers 树.R",
          paste0('XGBoost-tree/','0 Model RData/','XGBoost筛选biomarkers 树.R'),
          overwrite = T)









