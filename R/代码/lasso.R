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
      legend.text = element_text(face = "bold", color = "black", size = 20,family = 'Times'),
      legend.title = element_text(face = "bold", color = "black", size = 24,family = 'Times'),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'black', size = 1),
      axis.ticks = element_line(size = 1.2,linetype=2),
      axis.ticks.length = unit(0.3,'cm'),
      axis.text.y= element_text(face = "bold", color = data$col, size = 20),
      axis.text.x= element_text(face = "bold", color = 'black', size = 20),
      axis.title= element_text(face = "bold", color = "black", size = 24),
      title = element_text(face = "bold", color = "black", size = 24)
    )
  ggsave(
    paste0(output, 'varImpPlot.jpg'), p, 
    width = 10, height = max(0.2*nrow(data),8), units = 'in', dpi = 300
  )
  ggsave(
    paste0(output, 'varImpPlot.pdf'), p, 
    width = 10, height = max(0.2*nrow(data),8), units = 'in', dpi = 300
  )
}

rt=read.xlsx('./Random-forest/biomarker.xlsx',rowNames = T)
rt$group<-factor(rt$group,levels = unique(rt$group))
######一个简单的lasso
x=as.matrix(rt[,1:11])
y=gsub("(.*)\\_(.*)", "\\2", rt[,12])
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
#绘制Lasso回归的图形
pdf(file="lasso.pdf", width=6, height=5.5)
plot(fit, xvar = "lambda")
dev.off()
#绘制交叉验证的图形
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#输出筛选的特征基因
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]

gene_coef <- data.frame(biomarkers=lassoGene,coefficience=coef[index][-1])
write.csv(gene_coef,'./coefficience.csv')
importanceplot(gene_coef,output = './')
save(cvfit,file='Lasso.RData')


predata<-data.frame(y=predict(cvfit,newx = as.matrix(rt[,-ncol(rt)]),type='response'),group=rt$group)
colnames(predata)[1]<-'prediction'

predata<-predata[order(predata$prediction),]

predata$sample<-seq(0,1,length.out = nrow(rt))
p <- ggplot(data = predata, mapping = aes(x = sample,y = prediction,color=group))+
  geom_point(size=9)+
  theme_bw()+
  theme(
    legend.text = element_text(face = "bold", color = "black", size = nrow(rt)*0.36),
    legend.title = element_text(face = "bold", color = "black", size = nrow(rt)*0.56),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'black', size = nrow(rt)*0.03),
    axis.ticks = element_blank(),
    axis.text= element_blank(),
    axis.title = element_text(face = "bold", color = "black", size = nrow(rt)*0.56)
  )
ggsave(file='logistic regression.pdf',p,width = 30,height = 30)
ggsave(file='logistic regression.png',p,width = 30,height = 30,dpi=600)


rocdata<-data.frame(predict(cvfit,newx = as.matrix(rt[,-ncol(rt)]),type='class'))
colnames(rocdata)<-'prediction'
rocdata$true<-rt$group
write.csv(rocdata,file='logistic regression.csv')


if (length(lassoGene)>1) {
  set.seed(1234)
  ###用biomarkers重新绘制新模型
  cv.fittest<-cv.glmnet(as.matrix(rt[,lassoGene]),rt$group,family="binomial",type.measure='class',nfolds=10,standardize=F)
  save(cv.fittest,file='Biomarkers.RData')
  gt.coeftest <- coef(cv.fittest$glmnet.fit, s = cv.fittest$lambda.1se)
  
  #ROC曲线绘制
  if (!dir.exists(paste0('Lasso-logistic/','4 ROC curve'))){dir.create(paste0('Lasso-logistic/','4 ROC curve'))}
  getclassification<-function(data,output){
    classification<-data.frame(y=predict(cv.fittest,newx = as.matrix(data[,lassoGene]),type='class'))
    colnames(classification)[1]<-'prediction'
    classification$true<-data$group
    classification$prob<-unlist(predict(cv.fittest,newx = as.matrix(data[,lassoGene]),type='response'))
    write.csv(classification,file=output)
  }
  # #训练集
  rtrocs<-roc(as.numeric(rt$group)-1,predict(cv.fittest,newx = as.matrix(rt[,lassoGene]),type='response'),ci=T,percent=TRUE)
  getclassification(rt,output=paste0('Lasso-logistic/','4 ROC curve/','rt.csv'))
  # #验证集
  testdatarocs<-roc(as.numeric(testdata$group)-1,predict(cv.fittest,newx = as.matrix(testdata[,lassoGene]),type='response'),ci=T,percent=TRUE)
  getclassification(testdata,output=paste0('Lasso-logistic/','4 ROC curve/','testdata.csv'))
  # #得到新样本，观察模型预测新样本能力
  # 女性离群样本
  # liqundataraw<-read.xlsx('data_yz.xlsx',rowNames = T)
  # liqundataraw[,colnames(rt)[-ncol(rt)]]<-predict(preProc,liqundataraw[,colnames(rt)[-ncol(rt)]])#用训练集的均值和标准差来对新样本进行标准化
  # liqundata<-liqundataraw[, lassoGene]
  # liqundata$group<- unlist(sapply(rownames(liqundata),function(i){strsplit(i,'_')[[1]][1]}))
  # liqundata$group<- unlist(sapply(liqundata$group,function(i){ifelse(i=='C','Control',i)}))
  # liqundata$group<-factor(liqundata$group,levels = unique(liqundata$group))
  # 
  # if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
  #   liqundata$group<-factor(liqundata$group,levels = unique(liqundata$group))
  #   liqundataroc<-roc(as.numeric(liqundata$group)-1,predict(cv.fittest,newx = as.matrix(liqundata[,lassoGene]),type='response'),ci=T,percent=TRUE)
}

png(paste0('Lasso-logistic/','4 ROC curve/','ROC plot.png'))
plot.roc(rtrocs,print.thres=F,print.auc=F,col='red',reuse.auc=F)
# plot.roc(testdatarocs,print.thres=F,print.auc=F,col='blue',reuse.auc=F,add = T)
legendtest<-c('DATA: AUC(95%CI) %',
              paste0('rt: ',sprintf("%.2f (%.2f-%.2f) %%", rtrocs$auc, rtrocs$ci[1],rtrocs$ci[3]))
              # paste0('testdata: ',sprintf("%.2f (%.2f-%.2f) %%", testdatarocs$auc, testdatarocs$ci[1],testdatarocs$ci[3]))
)
legendcol<-c('black','red')#,'blue'
# if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
#   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
#   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
#   legendcol<-c(legendcol,'green')
# }
legend('bottomright',
       legendtest,col=legendcol,text.col = legendcol,bty = 'n')
dev.off()

Risk_score <- function(gene_coef,data){
  riskdata <- data
  riskdata$risk_score <- 0
  for(i in gene_coef$biomarkers){
    riskdata$risk_score <- riskdata$risk_score+data[,i]*gene_coef[gene_coef$biomarkers==i,"coefficience"]
  }
  riskdata <- riskdata[,c(match(gene_coef$biomarkers,colnames(riskdata)),ncol(riskdata),ncol(riskdata)-1)]
  return(riskdata)
}

filelist <- c('./01.差异分析/GSE32280.csv','./01.差异分析/GSE76826.csv','./01.差异分析/GSE98793.csv') 
for(i in filelist){
  data <- read.csv(i,row.names = 1)
  data <- as.data.frame(t(data))
  for(ii in colnames(data)[-ncol(data)]){
    data[,ii] <- as.numeric(data[,ii])
  }
  riskdata_ <- Risk_score(gene_coef = gene_coef,data = data)
  write.csv(riskdata_,paste0('risk_',strsplit(i,'/')[[1]][3]))
}

