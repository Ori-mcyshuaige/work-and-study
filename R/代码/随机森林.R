library(randomForest)
library(stringr)
library(ggpubr)
library(pals)###jet颜色配色
#library(doParallel) # 并行处理的包
library(Cairo)
library(openxlsx)
library(pROC)
#library(MUVR) # MUVR包

getvarImpPlot<-function(data,xlab,nvar=NULL,output=''){
  options(stringsAsFactors = F)

  getcolorname_scale_color_gradient<-function(x,low='blue',high='red', space = "Lab"){
    require(scales)
    get_palette<-seq_gradient_pal(low,high, space)
    get_rescaler<-function(x=x,to=c(0,1),from=range(x,na.rm=TRUE)){rescale(x,to,from)}
    return(get_palette(get_rescaler(x)))
  }
  
  data$Score<-as.numeric(data[,xlab])
  data<-data[order(data$Score,decreasing = F),]
  if(!is.null(nvar)){data<-data[((nrow(data)-nvar)+1):nrow(data),]}
  data$id<-factor(rownames(data),levels = rownames(data))
  data$col<-getcolorname_scale_color_gradient(data$Score)
  
  #show_col(data$col)
  p<- ggplot(data,aes(x=Score,y=id))+
    geom_segment(aes(yend=id,color=Score),xend=0,size=0.7)+
    geom_point(aes(color=Score),size=3)+
    #facet_grid(type~ .,scales = 'x',space='x')+
    scale_color_gradient(low='blue',high='red',name = xlab)+
    theme_bw() +
    theme(
      #legend.position = 'right', 
      #legend.justification = c(1,0),
      legend.text = element_text(face = "bold", color = "black", size = 20,family = 'Times'),
      legend.title = element_text(face = "bold", color = "black", size = 24,family = 'Times'),
      #legend.key.height = unit(5,'cm'),
      
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'black', size = 1),
      # axis.line.x.bottom=element_line(color = 'black', size = 1),
      # axis.line.y.left=element_line(color = 'black', size = 1),
      axis.ticks = element_line(size = 1.2,linetype=2),
      axis.ticks.length = unit(0.3,'cm'),
      axis.text.y= element_text(face = "bold", color = data$col, size = 12),
      axis.text.x= element_text(face = "bold", color = 'black', size = 12),
      axis.title= element_text(face = "bold", color = "black", size = 24)
    ) +
    labs(x = xlab)
  ggsave(
    paste0(output,xlab, '_varImpPlot.jpg'), p, 
    width = 12, height = max(0.02*nrow(data),8), units = 'in',limitsize = FALSE
  )
  ggsave(
    paste0(output,xlab, '_varImpPlot.pdf'), p, 
    width = 12, height = max(0.02*nrow(data),8), units = 'in',limitsize = FALSE
  )
}

scatterplot<-function(i,data,df.c){
  x <- df.c$x[i]
  y <- df.c$y[i]
  library(ggpubr)
  if (which(colnames(data)==x)<which(colnames(data)==y)) {
    p <- ggscatter(data, x = x,y = y,shape = 20,size = 5,
                   color = "group",
                   palette = jet(length(levels(data$group))+2)[2:(length(levels(data$group))+1)][as.numeric(unique(data$group))],add.params = list(color = "black"),
                   add = "reg.line",cor.method='spearman',
                   conf.int = T, cor.coef = T,show.legend.text = F)+
      geom_rug(aes(colour=group))+
      #scale_color_manual(values=jet(length(levels(data$group))+2)[2:(length(levels(data$group))+1)][as.numeric(unique(data$group))])+
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
    datax<-data[,c(x,"group")]
    datax$meta<-x
    colnames(datax)[1]<-'exep'
    datay<-data[,c(y,"group")]
    datay$meta<-y
    colnames(datay)[1]<-'exep'
    data<-rbind(datax,datay)
    #colnames(data)[1]<-'exep'
    ###均值
    descdata1<-as.data.frame(aggregate(data[,'exep'],list(meta=data[,'meta'],group=data[,'group']),mean))
    ##标准误
    descdata2<-as.data.frame(aggregate(data[,'exep'],list(meta=data[,'meta'],group=data[,'group']),function(i){sd(i)/sqrt(length(i))}))
    ##95%置信空间
    descdata3<-as.data.frame(aggregate(data[,'exep'],list(meta=data[,'meta'],group=data[,'group']),function(i){(sd(i)/sqrt(length(i)))*qt(.95/2 + .5, length(i)-1)}))
    descdata<-cbind.data.frame(descdata1,descdata2[,3],descdata3[,3])
    colnames(descdata)[3:5]<-c('mean','se','ci')
    data<-merge(data,descdata,by=c('meta','group'))
    
    p<-ggplot(descdata,aes(x=group,y=mean,fill=group,colour=meta,group=meta,linetype=meta))+
      geom_bar(stat = 'identity',position=position_dodge(0.4),width=.3)+
      #geom_point(data=data,aes(y=exep),position=position_dodge(0.4),size=1.2,alpha=0.5)+
      geom_line(position=position_dodge(0.4),size=.56)+
      geom_errorbar(aes(ymin=mean-ci,ymax=mean+ci),position=position_dodge(0.4),size=.56,width=.36)+
      annotate('text',label=unique(descdata$meta)[1],x=2,y=max(descdata$mean)*.95,size=12*.95,face='blod',color=jet(2)[1])+
      annotate('text',label=unique(descdata$meta)[2],x=2,y=max(descdata$mean)*.85,size=12*.95,face='blod',color=jet(2)[2])+
      scale_fill_manual(values = jet(length(levels(data$group))+2)[2:(length(levels(data$group))+1)][as.numeric(unique(data$group))])+
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
    descdata<-as.data.frame(aggregate(data[,'meta'],list(group=data[,'group']),mean))
    
    p<-ggplot(data,aes(x=as.numeric(group)-.3,y=meta,fill=group,color=group,group=group))+
      geom_violin()+
      geom_boxplot(width=.1,notch=T,fill='white')+
      geom_dotplot(aes(x=as.numeric(group)+.3,group=group),binaxis='y',stackdir = 'center',method='histodot')+
      geom_rug(sides='left')+
      annotate('text',label=y,x=2,y=max(data[,y])*.95,size=12*.95,face='blod')+
      scale_fill_manual(values=jet(length(levels(data$group))+2)[2:(length(levels(data$group))+1)][as.numeric(unique(data$group))])+
      scale_color_manual(values=jet(length(levels(data$group))+2)[2:(length(levels(data$group))+1)][as.numeric(unique(data$group))])+
      scale_x_continuous(breaks=1:nlevels(data$group),labels =levels(data$group))+
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

cvplot<-function(data,output=''){
  getcolorname_scale_color_gradient<-function(x,low='blue',high='red', space = "Lab"){
    require(scales)
    get_palette<-seq_gradient_pal(low,high, space)
    get_rescaler<-function(x=x,to=c(0,1),from=range(x,na.rm=TRUE)){rescale(x,to,from)}
    return(get_palette(get_rescaler(x)))
  }
  data<-data[order(data$variables),]
  data$fill<-getcolorname_scale_color_gradient(data$mean,high='red',low = 'grey')
  data$col<-getcolorname_scale_color_gradient(data$mean,high='blue',low = 'grey')
  p<-ggplot(data,aes(x=variables,y=mean,fill=mean,colour=se))+
    geom_bar(stat = 'identity')+
    geom_line(size=1.56)+
    geom_errorbar(aes(ymin=mean-ci,ymax=mean+ci),size=1.56)+
    scale_fill_gradient(low='grey',high='red')+
    scale_color_gradient(low='grey',high='blue')+
    labs(x='Number of Variables',y='CV Error')+
    annotate('text',x = max(data$variables)/2,y= max(data$mean)*0.95,size=9,
             label = 'Statistical Results of 10 Fold Cross Validation after 10 Times Repetition')+
    theme_bw()+
    theme(legend.position='right',
          legend.text = element_text(face = "bold", color = "black", size = 24,family = 'Times'),
          legend.title = element_text(face = "bold", color = "black", size = 32,family = 'Times'),
          
          axis.title = element_text(face = "bold",size = 40),
          #axis.ticks= element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks.length = unit(.56, "cm"),
          axis.text.x = element_text(face = "bold",color = data$fill,size = 36),
          axis.ticks.x=element_line(color = data$col,size=5),
          axis.text.y = element_text(face = "bold",size = 36),
          axis.ticks.y=element_line(size=5),
          legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = 'black', size = 1))
  ggsave(paste0(output,'CV-var Plot.png'),p,width = 20, height = 12, units = 'in', dpi = 300)
  ggsave(paste0(output,'CV-var Plot.pdf'),p,width = 20, height = 12, units = 'in', dpi = 300)
}





####导入数据
# raw<-read.xlsx('./01.差异分析/RF.xlsx',rowNames = T,sheet = 1,check.names = T,sep.names = " ")
# data<-raw
raw <- read.csv('./01.差异分析/df_diff_GSE98793.csv',row.names = 1)
data<-as.data.frame(t(raw))
for(i in colnames(data)[-ncol(data)]){
  data[,i] <- as.numeric(data[,i])
}
# data<-as.data.frame(t(raw))#[, -c(ncol(raw),ncol(raw)-1)]
# data$group<- unlist(sapply(rownames(data),function(i){strsplit(i,'_')[[1]][1]}))
# data$group<- unlist(sapply(data$group,function(i){ifelse(i=='C','Control',i)}))
data$group<-factor(data$group,levels = unique(data$group))
set.seed(16)
smin<-unlist(sapply(unique(data$group),function(i){sum(data$group==i)}))
names(smin)<-unique(data$group)
sid<-c(sample(rownames(data)[data$group==names(smin)[1]],3/4*smin[1]),sample(rownames(data)[data$group==names(smin)[2]],3/4*smin[2]))
sid<-unlist(sapply(sid, function(i){which(rownames(data)==i)}))
# sid<-sample(nrownames(data),3/4*nrow(data))
traindata<- data[sid,]
testdata<- data[-c(sid),]
for(ii in 1:ncol(traindata[,-ncol(traindata)])){
  result<-sapply(as.numeric(traindata[,ii]),function(i){(i-mean(as.numeric(traindata[,ii])))/sd(as.numeric(traindata[,ii]))})
  traindata[,ii]<-result
}
for(ii in 1:ncol(testdata[,-ncol(testdata)])){
  result<-sapply(as.numeric(testdata[,ii]),function(i){(i-mean(as.numeric(testdata[,ii])))/sd(as.numeric(testdata[,ii]))})
  testdata[,ii]<-result
}
# smin<-min(unlist(sapply(unique(train1data$group),function(i){sum(train1data$group==i)})))
# sid<-c(unlist(sapply(unique(train1data$group),function(i){sample(rownames(train1data)[train1data$group==i],3/4*smin)})))
# sid<-unlist(sapply(sid, function(i){which(rownames(train1data)==i)}))
# traindata<-train1data[sid,]
# validationdata<-train1data[-c(sid),]

# data1<-read.xlsx('RF.xlsx',rowNames = T,sheet = 2,check.names = T,sep.names = " ")
# data2<-read.xlsx('RF.xlsx',rowNames = T,sheet = 3,check.names = T,sep.names = " ")
# data3<-read.xlsx('RF.xlsx',rowNames = T,sheet = 4,check.names = T,sep.names = " ")
# traindata<-data[rownames(data1),]
# validationdata<-data[rownames(data2),]
# testdata<-data[rownames(data3),]

set.seed(16)
###随机森林模型
if (!dir.exists('Random-forest')){dir.create('Random-forest')}
data.rf = randomForest(group ~ ., data=traindata, importance=TRUE, proximity=TRUE,ntree=500)
if (!dir.exists(paste0('Random-forest/','0 Model RData'))){dir.create(paste0('Random-forest/','0 Model RData'))}
save(data.rf,file=paste0('Random-forest/','0 Model RData/','randomForest.RData'))
write.csv(traindata,paste0('Random-forest/','0 Model RData/','traindata.csv'),row.names =T)
# write.csv(validationdata,paste0('Random-forest/','0 Model RData/','validationdata.csv'),row.names=T)
write.csv(testdata,paste0('Random-forest/','0 Model RData/','testdata.csv'),row.names=T)


bioim=as.data.frame(data.rf$importance)
varim<-as.data.frame(importance(data.rf,class = NULL, scale = TRUE, type = NULL))
varim<-bioim
if (!dir.exists(paste0('Random-forest/','1 metabolite_importance'))){dir.create(paste0('Random-forest/','1 metabolite_importance'))}
write.csv(varim,paste0('Random-forest/','1 metabolite_importance/','metabolite_importance.csv'))
###绘图varImpPlot重新优化 新函数getvarImpPlot
if (!dir.exists(paste0('Random-forest/','1 metabolite_importance/all'))){dir.create(paste0('Random-forest/','1 metabolite_importance/all'))}
getvarImpPlot(varim,'MeanDecreaseAccuracy',output=paste0('Random-forest/','1 metabolite_importance/','all/'))
getvarImpPlot(varim,'MeanDecreaseGini',output=paste0('Random-forest/','1 metabolite_importance/','all/'))



###分类图MDSplot 待优化
#rf:随机森林对象; fac:训练随机森林的因子; k:维度数; palette:颜色
if (!dir.exists(paste0('Random-forest/','2 MDS Plot'))){dir.create(paste0('Random-forest/','2 MDS Plot'))}
png(file=paste0("MDS Plot.png"),width=560,height=560)
MDSplot(rf=data.rf,fac=traindata$group,k=2,#pch = as.numeric(data$group),
        palette=jet(nlevels(traindata$group)+2)[2:(nlevels(traindata$group)+1)][as.numeric(unique(traindata$group))])
dev.off()
rf_mds<-stats::cmdscale(1 - data.rf$proximity, eig = TRUE,k = 2)
colnames(rf_mds$points) <- paste("Dim", 1:2)
write.csv(rf_mds$points,paste0('Random-forest/','2 MDS Plot/',"MDS Plot.csv"),row.names = T)
pdf(file=paste0("MDS Plot.pdf"),width=5.6,height=5.6)
MDSplot(rf=data.rf,fac=traindata$group,k=2,#pch = as.numeric(data$group),
        palette=jet(nlevels(traindata$group)+2)[2:(nlevels(traindata$group)+1)][as.numeric(unique(traindata$group))])
dev.off()
## 无监督分类
###判断各组间是否存在差异
set.seed(16)
data.urf= randomForest(traindata[, -ncol(traindata)],importance=TRUE, proximity=TRUE,ntree=500)
# 主坐标轴分析并展示
png(file=paste0("Unsupervised MDS Plot.png"),width=560,height=560)
MDSplot(rf=data.urf,fac=traindata$group,k=2,palette=jet(nlevels(traindata$group)+2)[2:(nlevels(traindata$group)+1)][as.numeric(unique(traindata$group))])
dev.off()
pdf(file=paste0("Unsupervised MDS Plot.pdf"),width=5.6,height=5.6)
MDSplot(rf=data.urf,fac=traindata$group,k=2,palette=jet(nlevels(traindata$group)+2)[2:(nlevels(traindata$group)+1)][as.numeric(unique(traindata$group))])
dev.off()
urf_mds<-stats::cmdscale(1 - data.urf$proximity, eig = TRUE,k = 2)
colnames(urf_mds$points) <- paste("Dim", 1:2)
write.csv(urf_mds$points,paste0('Random-forest/','2 MDS Plot/',"Unsupervised MDS Plot.csv"),row.names = T)

for (plotfile in c("MDS Plot.png","MDS Plot.pdf","Unsupervised MDS Plot.png","Unsupervised MDS Plot.pdf")) {
  file.copy(plotfile,paste0('Random-forest/','2 MDS Plot/',plotfile),overwrite=T)
  file.remove(plotfile) 
}




##交叉验证进行特征选择：rfcv()
if (!dir.exists(paste0('Random-forest/','3 CVerror'))){dir.create(paste0('Random-forest/','3 CVerror'))}
set.seed(16)
result <- replicate(10, rfcv(traindata[,-ncol(traindata)], traindata$group, cv.fold=10), simplify=FALSE)
error.cv <- sapply(result, "[[", "error.cv")
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l",
        lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")

cvdata<-data.frame(mean=sapply(1:nrow(error.cv), function(i){mean(error.cv[i,])}),
                   se=sapply(1:nrow(error.cv), function(i){sd(error.cv[i,])/sqrt(ncol(error.cv))}),
                   ci=sapply(1:nrow(error.cv), function(i){(sd(error.cv[i,])/sqrt(ncol(error.cv)))*qt(.95/2 + .5, ncol(error.cv)-1)})
)
rownames(cvdata)<-rownames(error.cv)
cvdata$variables<-as.numeric(rownames(cvdata))
cvplot(cvdata,output = paste0('Random-forest/','3 CVerror/'))
write.csv(cbind.data.frame(cvdata,error.cv),file=paste0('Random-forest/','3 CVerror/CVerror.csv'),row.names = F)

###选择最优变量个数
repeat{
  nvar<-cvdata$variables[which.min(cvdata$mean)]
  if (nvar<=10) {
    break()
  }else{
    cat(getwd())
    nvar<- as.numeric(readline('Number of variables: '))
  }
  if(nvar <= 30){break()}
  else{
    cat(getwd())
    ff<- readline('Number of variables is larger. Continue?: Y or N  : ')
    if (ff=='Y') {
      break()
      }
  }
}

###选取最优变量组合
# bios=rownames(varim[order(varim$MeanDecreaseAccuracy,decreasing = T)[1:nvar],])
bios=intersect(rownames(varim[order(varim$MeanDecreaseAccuracy,decreasing = T)[1:nvar],]),
               rownames(varim[order(varim$MeanDecreaseGini,decreasing = T)[1:nvar],]))
gooddata1<-traindata[,c(bios,'group')]
# gooddata2<-validationdata[,c(bios,'group')]
# gooddata3<-testdata[,c(bios,'group')]
# write.xlsx(gooddata2,paste0('Random-forest/','biomarker-validation.xlsx'),rowNames=T)
# write.xlsx(gooddata1,paste0('Random-forest/','biomarker-train.xlsx'),rowNames=T)
# write.xlsx(gooddata3,paste0('Random-forest/','biomarker-test.xlsx'),rowNames=T)
gooddata<-data[,c(bios,'group')]
write.xlsx(gooddata,paste0('Random-forest/','biomarker.xlsx'),rowNames=T)
#goodvarim<-varim[bios,]
###绘图varImpPlot重新优化 新函数getvarImpPlot
if (!dir.exists(paste0('Random-forest/','1 metabolite_importance/biomarker'))){dir.create(paste0('Random-forest/','1 metabolite_importance/biomarker'))}
getvarImpPlot(varim,'MeanDecreaseAccuracy',nvar = nvar,output=paste0('Random-forest/','1 metabolite_importance','/biomarker/'))
getvarImpPlot(varim,'MeanDecreaseGini',nvar = nvar,output=paste0('Random-forest/','1 metabolite_importance','/biomarker/'))


###画变量关系图
if (!dir.exists(paste0('Random-forest/','4 Features plot'))){dir.create(paste0('Random-forest/','4 Features plot'))}
p<-Features(gooddata,output = paste0('Random-forest/','4 Features plot/'))
CairoPNG(file="Features plot.png",width=(ncol(traindata)+1)*56,height=(ncol(traindata)+1)*56)
grid.draw(p)
dev.off()
file.copy("Features plot.png",paste0('Random-forest/','4 Features plot/','Features plot.png'),overwrite=T)
file.remove("Features plot.png")


###重构模型，验证集ROC

if (nvar>1) {
  if (!dir.exists(paste0('Random-forest/','5 ROC curve'))){dir.create(paste0('Random-forest/','5 ROC curve'))}
  set.seed(16)
  data.rfnew = randomForest(group ~ ., data=gooddata1, importance=TRUE, proximity=TRUE,ntree=500)
  save(data.rfnew,file=paste0('Random-forest/','0 Model RData/','Biomarkers.RData'))
  
  getclassification<-function(data,output){
    classification<-data.frame(y=predict(data.rfnew,data[,bios],type='class'))
    colnames(classification)[1]<-'prediction'
    classification$true<-data$group
    classification$prob<-predict(data.rfnew,data[,bios],type='prob')[,2]
    write.csv(classification,file=output)
  }
  # #训练集
  traindatarocs<-roc(as.numeric(traindata$group)-1,predict(data.rfnew,traindata[,bios],type='prob')[,2],ci=T,percent=TRUE)
  getclassification(traindata,output=paste0('Random-forest/','5 ROC curve/','traindata.csv'))
  # #验证集1
  # validationdatarocs<-roc(as.numeric(validationdata$group)-1,predict(data.rfnew,validationdata[,bios],type='prob')[,2],ci=T,percent=TRUE)
  # getclassification(validationdata,output=paste0('Random-forest/','5 ROC curve/','validationdata.csv'))
  # # #得到新样本，观察模型预测新样本能力
  # # #验证集2
  # liqundataraw<-read.xlsx('data_yz.xlsx',rowNames = T)
  # liqundata<-liqundataraw[, bios]
  # # liqundata<-data.frame(scale(liqundata))
  # liqundata$group<- unlist(sapply(rownames(liqundata),function(i){strsplit(i,'_')[[1]][1]}))
  # liqundata$group<- unlist(sapply(liqundata$group,function(i){ifelse(i=='C','Control',i)}))
  # liqundata$group<-factor(liqundata$group,levels = unique(liqundata$group))
  # 
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   liqundata$group<-factor(liqundata$group,levels = unique(liqundata$group))
  #   liqundataroc<-roc(as.numeric(liqundata$group)-1,predict(data.rfnew,liqundata[,bios],type='prob')[,2],ci=T,percent=TRUE)
  # }
  # getclassification(liqundata,output=paste0('Random-forest/','5 ROC curve/','离群.csv'))
  testdatarocs<-roc(as.numeric(testdata$group)-1,predict(data.rfnew,testdata[,bios],type='prob')[,2],ci=T,percent=TRUE)
  getclassification(testdata,output=paste0('Random-forest/','5 ROC curve/','testdata.csv'))
  
  
  png(paste0('Random-forest/','5 ROC curve/','ROC plot.png'),width = 400,height = 400)
  plot.roc(traindatarocs,ci=T,plot=T,grid=F,print.auc=T,legacy.axes=T,print.thres="best",thresholds="best",reuse.auc=FALSE)
  # plot.roc(validationdatarocs,print.thres=F,print.auc=F,col='#FFD400',reuse.auc=F,add = T)
  plot.roc(testdatarocs,ci=T,plot=T,grid=F,print.auc=T,legacy.axes=T,print.thres="best",thresholds="best",reuse.auc=FALSE,add = T)
  legendtest<-c('DATA: AUC(95%CI) %',
                paste0('train: ',sprintf("%.2f (%.2f-%.2f) %%", traindatarocs$auc, traindatarocs$ci[1],traindatarocs$ci[3])),
                # paste0('test1: ',sprintf("%.2f (%.2f-%.2f) %%", validationdatarocs$auc, validationdatarocs$ci[1],validationdatarocs$ci[3]))
                paste0('test: ',sprintf("%.2f (%.2f-%.2f) %%", testdatarocs$auc, testdatarocs$ci[1],testdatarocs$ci[3]))
  )
  legendcol<-c("black",'#00D4FF','#7F0000')
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
  #   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
  #   legendcol<-c(legendcol,'green')
  # }
  legend('bottomright',
         legendtest,col=legendcol,text.col = legendcol,bty = 'n')
  
  
  dev.off()
  
  pdf(paste0('Random-forest/','5 ROC curve/','ROC plot.pdf'),width = 5.5,height = 5.5)
  plot.roc(traindatarocs,ci=T,plot=T,grid=F,print.auc=T,legacy.axes=T,print.thres="best",thresholds="best",reuse.auc=FALSE)
  # plot.roc(validationdatarocs,print.thres=F,print.auc=F,col='#FFD400',reuse.auc=F,add = T)
  plot.roc(testdatarocs,ci=T,plot=T,grid=F,print.auc=T,legacy.axes=T,print.thres="best",thresholds="best",reuse.auc=FALSE,add = T)
  legendtest<-c('DATA: AUC(95%CI) %',
                paste0('train: ',sprintf("%.2f (%.2f-%.2f) %%", traindatarocs$auc, traindatarocs$ci[1],traindatarocs$ci[3])),
                # paste0('test1: ',sprintf("%.2f (%.2f-%.2f) %%", validationdatarocs$auc, validationdatarocs$ci[1],validationdatarocs$ci[3]))
                paste0('test: ',sprintf("%.2f (%.2f-%.2f) %%", testdatarocs$auc, testdatarocs$ci[1],testdatarocs$ci[3]))
  )
  legendcol<-c("black",'#00D4FF','#7F0000')
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
  #   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
  #   legendcol<-c(legendcol,'green')
  # }
  legend('bottomright',
         legendtest,col=legendcol,text.col = legendcol,bty = 'n')
  
  
  dev.off()
  
  
  pdf(paste0('Random-forest/','5 ROC curve/','train ROC plot.pdf'),width = 5.5,height = 5.5)
  plot.roc(traindatarocs,ci=T,plot=T,grid=F,print.auc=T,legacy.axes=T,print.thres="best",thresholds="best",reuse.auc=FALSE)
  legendtest<-c('DATA: AUC(95%CI) %',
                paste0('train: ',sprintf("%.2f (%.2f-%.2f) %%", traindatarocs$auc, traindatarocs$ci[1],traindatarocs$ci[3]))
  )
  legendcol<-c('black','black')
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
  #   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
  #   legendcol<-c(legendcol,'green')
  # }
  legend('bottomright',
         legendtest,col=legendcol,text.col = legendcol,bty = 'n')
  
  
  dev.off()
  
  # pdf(paste0('Random-forest/','5 ROC curve/','test1 ROC plot.pdf'),width = 5.5,height = 5.5)
  # plot.roc(validationdatarocs,print.thres=F,print.auc=F,col='black',reuse.auc=F)
  # legendtest<-c('DATA: AUC(95%CI) %',
  #               paste0('test1: ',sprintf("%.2f (%.2f-%.2f) %%", validationdatarocs$auc, validationdatarocs$ci[1],validationdatarocs$ci[3]))
  # )
  # legendcol<-c('black','black')
  # # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  # #   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
  # #   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
  # #   legendcol<-c(legendcol,'green')
  # # }
  # legend('bottomright',
  #        legendtest,col=legendcol,text.col = legendcol,bty = 'n')
  # 
  # 
  # dev.off()
  
  pdf(paste0('Random-forest/','5 ROC curve/','test ROC plot.pdf'),width = 5.5,height = 5.5)
  plot.roc(testdatarocs,ci=T,plot=T,grid=F,print.auc=T,legacy.axes=T,print.thres="best",thresholds="best",reuse.auc=FALSE)
  legendtest<-c('DATA: AUC(95%CI) %',
                paste0('test2: ',sprintf("%.2f (%.2f-%.2f) %%", testdatarocs$auc, testdatarocs$ci[1],testdatarocs$ci[3]))
  )
  legendcol<-c('black','black')
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
  #   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
  #   legendcol<-c(legendcol,'green')
  # }
  legend('bottomright',
         legendtest,col=legendcol,text.col = legendcol,bty = 'n')
  
  
  dev.off()
  save(traindatarocs,file=paste0('Random-forest/','5 ROC curve/','trainrocPRO.RData'))
  # save(validationdatarocs,file=paste0('Random-forest/','5 ROC curve/','validationrocmet.RData'))
  save(testdatarocs,file=paste0('Random-forest/','5 ROC curve/','testrocPRO.RData'))

  
  # rocs<-roc(as.numeric(testdata$group)-1,as.numeric(predict(data.rfnew,testdata[,tzid]))-1,ci=T,percent=TRUE)
  # png(paste0('Random-forest/','5 ROC curve/','ROC plot.png'))
  # plot.roc(rocs,print.thres=TRUE,print.auc=TRUE,col='red',reuse.auc=FALSE)
  # dev.off()
  # pdf(paste0('Random-forest/','5 ROC curve/','ROC plot.pdf'))
  # plot.roc(rocs,print.thres=TRUE,print.auc=TRUE,col='red',reuse.auc=FALSE)
  # dev.off()
  # 
  # 
  # ###新的样本验证模型
  # validationraw<-read.csv('data_yz.csv',row.names = 1)
  # if (length(unique(unlist(sapply(rownames(validationraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))==1){
  #   set.seed(16)
  #   validationraw<-rbind(validationraw, raw[sample(rownames(raw)[data$group=='HB'],sum(data$group=='HB')/2),])
  # }
  # validationdata<-validationraw[, tzid]
  # # validationdata<-data.frame(log(validationdata))
  # validationdata$group<- unlist(sapply(rownames(validationdata),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))
  
  validationdata$group<-factor(validationdata$group,levels = unique(validationdata$group))
  
  rocs<-roc(as.numeric(validationdata$group)-1,predict(data.rfnew,validationdata[,-ncol(validationdata)],type='prob')[,2],ci=T,percent=TRUE)
  getclassification(validationdata,output=paste0('Random-forest/','5 ROC curve/','validationdata.csv'))
  png('validationdata samples ROC plot.png')
  plot.roc(rocs,print.auc=TRUE,col='red',reuse.auc=FALSE,plot=T,grid=F,legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
  pdf('validationdata samples ROC plot.pdf')
  plot.roc(rocs,print.auc=TRUE,col='red',reuse.auc=FALSE,plot=T,grid=F,legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
}

validationdataall<-read.xlsx("Supplementary table S6-Targeted analysis.xlsx",sheet = 1,rowNames = T)
colnames(validationdataall)<-unlist(validationdataall[1,])
validationdata<-validationdataall[-1,c(2,3,6)]
for(ii in 1:(ncol(validationdata)-1)){
  result<-sapply(as.numeric(validationdata[,ii]),function(i){(i-mean(as.numeric(validationdata[,ii])))/sd(as.numeric(validationdata[,ii]))})
  validationdata[,ii]<-result
}
validationdata$group<-str_remove_all(rownames(validationdata),"\\d")
# 
# 
# validationdata<-validationdata[c(rocs$controls==0|1,rocs$cases==1),]
# ###另一个脚本
# file.copy("随机森林.R",paste0('Random-forest/','0 Model RData/','随机森林.R'),overwrite = T)
# #file.copy("biomarker.xlsx",paste0('Random-forest/','5 ROC curve/','biomarker.xlsx'),overwrite = T)
# 
# library(ggplot2)
# library(reshape2) #融合数据
# 
# #训练集混淆矩阵
# infile<-data.rf$confusion
# write.csv(infile,file = '训练集混淆矩阵表.csv',row.names = T)
# infile<-infile[,1:(ncol(infile)-1)]
# heatmapIN<-melt(as.matrix(infile))
# 
# 
# p<-ggplot(data = heatmapIN,
#           aes(x=Var1,y=Var2,fill=value)
# ) + geom_tile()+
#   scale_fill_gradient2(low = 'white',high = 'dodgerblue4',
#                        name=' '
#   )+theme_minimal()+
#   theme(axis.text.x = element_text(angle = 0,size = 15),axis.text.y = element_text(size = 15))+
#   labs(x="Predicted",y='True',title = "Train data")+geom_text(aes(Var1, Var2, label = value), color = "gray1", size = 5)+
#   theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+theme(plot.title = element_text(hjust = 0.5,size = 25))+
#   theme(legend.text = element_text(size = 12))+theme(legend.key.width=unit(3,'mm'),legend.key.height=unit(3,'cm'))
# 
# ggsave(p,filename = '训练集混淆矩阵热图.png',dpi=600,height = 9.6,width = 9.6)
# 
# 
# #验证集混淆矩阵
# infile<-data.rfnew$confusion
# write.csv(infile,file = '验证集混淆矩阵表.csv',row.names = T)
# infile<-infile[,1:(ncol(infile)-1)]
# heatmapIN<-melt(as.matrix(infile))
# 
# 
# p<-ggplot(data = heatmapIN,
#           aes(x=Var1,y=Var2,fill=value)
# ) + geom_tile()+
#   scale_fill_gradient2(low = 'white',high = 'dodgerblue4',
#                        name=' '
#   )+theme_minimal()+
#   theme(axis.text.x = element_text(angle = 0,size = 15),axis.text.y = element_text(size = 15))+
#   labs(x="Predicted",y='True',title = "Test data")+geom_text(aes(Var1, Var2, label = value), color = "gray1", size = 5)+
#   theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+theme(plot.title = element_text(hjust = 0.5,size = 25))+
#   theme(legend.text = element_text(size = 12))+theme(legend.key.width=unit(3,'mm'),legend.key.height=unit(3,'cm'))
# 
# ggsave(p,filename = '验证集混淆矩阵热图.png',dpi=600,height = 9.6,width = 9.6)

