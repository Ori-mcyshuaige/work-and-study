library(openxlsx)
library(stringr)
library(psych)
library(corrplot)
library(pheatmap)
library(ggpubr)
library(MASS)
library(readxl)


setwd("D:/马驰宇/桌面/何院专辑/研究院项目/戴主任项目重分析/关联分析/岭回归")


# corrall<-as.data.frame(matrix(nrow = 1,ncol = 1))
# pvalueall<-as.data.frame(matrix(nrow = 1,ncol = 1))
# for(sheet in excel_sheets("TGroup.xlsx")){
# temp0<-sheet
temp0<-"26代谢物重新分组"
if(!dir.exists(temp0)){dir.create(temp0)}

data1<-read.xlsx("data.xlsx",sheet = 2,colNames = T,rowNames = T,check.names = F,sep.names = ".")
data2<-read.xlsx("data.xlsx",sheet = 1,colNames = T,rowNames = T,check.names = F,sep.names = ".")
data2<-data2[rownames(data1),]
# pol<-read.xlsx("TGroup.xlsx",sheet = sheet,check.names = F,sep.names = " ")
# pos<-pol$POS[!is.na(pol$POS)]
# data2<-data2[rownames(data1),]
# data2<-data2[pos,]
# data1<-as.data.frame(data1[pos,])
# rownames(data1)<-pos
# colnames(data1)<-"Cholesterol"
# rownames(data1)<-str_replace(rownames(data1),"SLE",sheet)
# rownames(data2)<-str_replace(rownames(data2),"SLE",sheet)
corrdata<-corr.test(data2,data1[,1],method = 'spearman',adjust = 'none',use = 'p')
# datafilt<-rownames(corrdata$p)[corrdata$p<0.05]
datafilt<-rownames(corrdata$p)[corrdata$p<1]
data<-data2[,datafilt]
plotdata<-as.data.frame(t(data))
write.xlsx(plotdata,paste0(temp0,'/expression.xlsx'),rowNames=T)
class<-read.xlsx("data.xlsx",sheet = 3,colNames = T,rowNames = T,check.names = F,sep.names = ".")
class<-as.data.frame(t(class))
class<-as.data.frame(class[rownames(plotdata),])
# colnames(class)<-"group"
# rownames(class)<-rownames(plotdata)
corr<-data.frame(all=rep(1,length(colnames(data1))),row.names = colnames(data1))
pvalue<-data.frame(all=rep(1,length(colnames(data1))),row.names = colnames(data1))
  
for(i in colnames(data1)){
  data<-data2[,datafilt]
  data$result<-data1[,i]
  ridge.dia<-lm.ridge(result~.,lambda = seq(0,1000,length.out = 1001),data=data[!is.na(data$result),],model = TRUE)  #拟合回归模型
  # names(ridge.dia)  #查看回归模型的详细内容
  # "coef"   "scales" "Inter"  "lambda" "ym"     "xm"     "GCV"    "kHKB"   "kLW"
  # #coef：回归系数矩阵，每一行对应一个λ。
  # #scales：自变量矩阵的标度。
  # #Inter：是否包括截距？
  # #lambda：λ向量。
  # #ym：因变量的均值。
  # #xm：自变量矩阵按列计算的均值。
  # #GCV：广义交叉验证GVC向量。
  # #kHKB：岭常量的kHKB估计。
  # #kLW：岭常量的L-W估计。
  ridge.dia$lambda[which.min(ridge.dia$GCV)]	#输出使GCV最小时的λ值
  ridge.dia$coef[,which.min(ridge.dia$GCV)] #找到GCV最小时对应的系数，该结果是一个向量，长度为自变量的个数
  par(mfrow=c(1,2))
  #绘制回归系数关于λ的图形，并作出GCV取取最小值的那条竖线。
  matplot(ridge.dia$lambda,t(ridge.dia$coef),xlab = expression(lambda),ylab = "Coefficients",type = "l",lty = 1:20)
  abline(v=ridge.dia$lambda[which.min(ridge.dia$GCV)])
  #绘制GCV关于λ的图形，并作出GCV取取最小值的那条竖线。
  plot(ridge.dia$lambda,ridge.dia$GCV,type = "l",xlab = expression(lambda),ylab = expression(beta))
  abline(v=ridge.dia$lambda[which.min(ridge.dia$GCV)])
  
  mouse.regression<-lm.ridge(result~.,lambda = ridge.dia$lambda[which.min(ridge.dia$GCV)],data=data[!is.na(data$result),],model = TRUE)  #拟合回归模型
  
  
  
  
  
  con <- file(paste(temp0,'/',i,".log", sep = "")) # 创建一个.log文件
  sink(con, append=TRUE) # 记录output
  sink(con, append=TRUE, type="message") # 记录message
  # 所有的output和message都会记录到test.log中，而控制台中不在有信息显示
  
  # mouse.regression
  source('qqw.R')
  
  # 记录完毕后，重置output和message的记录，运行完一下两行，后续的输入命令重新显示到控制台中
  sink()
  sink(type="message")
  
  outprint<-readLines(paste(temp0,'/',i,".log", sep = ""))
  ns<-c('');xs<-c()
  for (ss in 1:length(outprint)) {
    sssss<-strsplit(outprint[ss],split = ' ')[[1]]
    sssss<-sssss[sssss!='']
    if(ss%%2==1){
      ns<-c(ns,sssss)
    }else{
      xs<-c(xs,sssss)
    }
  }
  xs<-sapply(xs, as.numeric)
  names(xs)<-sapply(ns, function(i){str_remove_all(i,'`')})
  names(xs)[1]<-'as'
  
  xs<-xs[c(colnames(data)[-ncol(data)],'as')]
  data$predict<-sapply(1:nrow(data),function(i){
    sum(as.numeric(data[i,-ncol(data)])*xs[-ncol(data)])+xs[ncol(data)]
  })
  
  p<-cor.test(data1[,i],data[rownames(data1),'predict'],method = 'spearman')$p.value
  c<-cor.test(data1[,i],data[rownames(data1),'predict'],method = 'spearman')$estimate
  corr[i,"all"]<-c
  pvalue[i,"all"]<-p
  
  data$Group<-unlist(sapply(rownames(data),function(i){paste0(strsplit(i,'_')[[1]][1],collapse = '')}))
  p<-ggscatter(data[!is.na(data$result),],'predict','result',cor.method ='spearman',shape = 20,size = 5,add = "reg.line",
               color = "Group",palette = c("#CD6090", "#1E90FF"),add.params = list(color = "black"),
               conf.int = T, cor.coef = T,show.legend.text = TRUE)+
    # geom_rug(aes(colour=Group))+
    theme_bw()+
    ylab('Result')+
    theme(
      #legend.position=c(.925,.5),
      #修改字号
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      axis.title = element_text(size = 20 ,colour = "black"),
      #end
      legend.key = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = 'black', size = 1)
    )
  
  ggsave(paste(temp0,'/',i,".png", sep = ""), p, width = 7.5, height = 7.5, units = 'in', dpi = 300)
  ggsave(paste(temp0,'/',i,".pdf", sep = ""),p,width = 7.5,height = 7.5,units = 'in',dpi = 300)
}


data2class<-class$group
for (jj in unique(data2class)) {
  temp1<-paste0(temp0,'/',jj)
  if(!dir.exists(temp1)){dir.create(temp1)}
  for (ii in colnames(data1)) {
    if(sum(data2class==jj)<2){
      data<-data2[,rep(datafilt[data2class==jj],3)]
      colnames(data)[c(2,3)]<-c('result','predict')
      data$result<-data1[,ii]
    }else{
      data<-data2[,datafilt[data2class==jj]]
      data$result<-data1[,ii]
      ## 岭回归
      ridge.dia<-lm.ridge(result~.,lambda = seq(0,1000,length.out = 1001),data=data[!is.na(data$result),],model = TRUE)  #拟合回归模型
      # names(ridge.dia)  #查看回归模型的详细内容
      # "coef"   "scales" "Inter"  "lambda" "ym"     "xm"     "GCV"    "kHKB"   "kLW"
      # #coef：回归系数矩阵，每一行对应一个λ。
      # #scales：自变量矩阵的标度。
      # #Inter：是否包括截距？
      # #lambda：λ向量。
      # #ym：因变量的均值。
      # #xm：自变量矩阵按列计算的均值。
      # #GCV：广义交叉验证GVC向量。
      # #kHKB：岭常量的kHKB估计。
      # #kLW：岭常量的L-W估计。
      ridge.dia$lambda[which.min(ridge.dia$GCV)]	#输出使GCV最小时的λ值
      ridge.dia$coef[,which.min(ridge.dia$GCV)] #找到GCV最小时对应的系数，该结果是一个向量，长度为自变量的个数
      par(mfrow=c(1,2))
      #绘制回归系数关于λ的图形，并作出GCV取取最小值的那条竖线。
      matplot(ridge.dia$lambda,t(ridge.dia$coef),xlab = expression(lambda),ylab = "Coefficients",type = "l",lty = 1:20)
      abline(v=ridge.dia$lambda[which.min(ridge.dia$GCV)])
      #绘制GCV关于λ的图形，并作出GCV取取最小值的那条竖线。
      plot(ridge.dia$lambda,ridge.dia$GCV,type = "l",xlab = expression(lambda),ylab = expression(beta))
      abline(v=ridge.dia$lambda[which.min(ridge.dia$GCV)])
      
      mouse.regression<-lm.ridge(result~.,lambda = ridge.dia$lambda[which.min(ridge.dia$GCV)],data=data[!is.na(data$result),],model = TRUE)  #拟合回归模型
      
      
      
      
      
      con <- file(paste(temp1,'/',ii,".log", sep = "")) # 创建一个.log文件
      sink(con, append=TRUE) # 记录output
      sink(con, append=TRUE, type="message") # 记录message
      # 所有的output和message都会记录到test.log中，而控制台中不在有信息显示
      
      source('qqw.R')
      
      # 记录完毕后，重置output和message的记录，运行完一下两行，后续的输入命令重新显示到控制台中
      sink()
      sink(type="message")
      
      outprint<-readLines(paste(temp1,'/',ii,".log", sep = ""))
      ns<-c('');xs<-c()
      for (ss in 1:length(outprint)) {
        sssss<-strsplit(outprint[ss],split = ' ')[[1]]
        sssss<-sssss[sssss!='']
        if(ss%%2==1){
          ns<-c(ns,sssss)
        }else{
          xs<-c(xs,sssss)
        }
      }
      xs<-sapply(xs, as.numeric)
      names(xs)<-sapply(ns, function(i){str_remove_all(i,'`')})
      names(xs)[1]<-'as'
      
      xs<-xs[c(colnames(data)[-ncol(data)],'as')]
      data$predict<-sapply(1:nrow(data),function(i){
        sum(as.numeric(data[i,-ncol(data)])*xs[-ncol(data)])+xs[ncol(data)]
      })
    }
      
      p<-cor.test(data1[,ii],data[rownames(data1),'predict'],method = 'spearman')$p.value
      c<-cor.test(data1[,ii],data[rownames(data1),'predict'],method = 'spearman')$estimate
      corr[ii,jj]<-c
      pvalue[ii,jj]<-p
        
    data$Group<-unlist(sapply(rownames(data),function(i){paste0(strsplit(i,'_')[[1]][1],collapse = '')}))
    p<-ggscatter(data[!is.na(data$result),],'predict','result',cor.method ='spearman',shape = 20,size = 5,add = "reg.line",
                 color = "Group",palette = c("#CD6090", "#1E90FF"),add.params = list(color = "black"),
                 conf.int = T, cor.coef = T,show.legend.text = TRUE)+
      # geom_rug(aes(colour=Group))+
      theme_bw()+
      ylab('Result')+
      theme(
        #legend.position=c(.925,.5),
        #修改字号
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 20 ,colour = "black"),
        #end
        legend.key = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 1)
      )
    
    ggsave(paste(temp1,'/',ii,".jpg", sep = ""), p, width = 6, height = 7.5, units = 'in', dpi = 600)
    ggsave(paste(temp1,'/',ii,".pdf", sep = ""),p,width = 6,height = 7.5,units = 'in',dpi = 600)
  }
}

# corr$group<-sheet
# pvalue$group<-sheet
# corrall<-rbind(corrall,corr)
# pvalueall<-rbind(pvalueall,pvalue)
# colnames(corrall)<-colnames(corr)
# colnames(pvalueall)<-colnames(pvalue)


# corrdata<-corr.test(data2,data1,method = 'spearman',adjust = 'none',use = 'p')
# 
# spearmanr<-data.frame(corrdata$r)
# spearmanp<-data.frame(corrdata$p)
# spearmanq
# # colnames(data)<-'corr'
# data$p_pasi<-sapply(rownames(data),function(i){cor.test(lipiddata[,i],clinical[,'Cholesterol'],method = 'spearman')$p.value})
# data$p_bas<-sapply(rownames(data),function(i){cor.test(lipiddata[,i],clinical[,'BAS'],method = 'spearman')$p.value})
# data$p_pga<-sapply(rownames(data),function(i){cor.test(lipiddata[,i],clinical[,'PGA'],method = 'spearman')$p.value})
# data$q_pasi<-p.adjust(data$p_pasi,method = 'BH')
# data$q_bas<-p.adjust(data$p_bas,method = 'BH')
# data$q_pga<-p.adjust(data$p_pga,method = 'BH')
# colnames(data)<-c('corr(PASI)','corr(BAS)','corr(PGA)','P-Value(PASI)','P-Value(BAS)','P-Value(PGA)','Q-Value(PASI)','Q-Value(BAS)','Q-Value(PGA)')
# 
# write.xlsx(data[,c('corr(PASI)','P-Value(PASI)','Q-Value(PASI)','corr(BAS)','P-Value(BAS)','Q-Value(BAS)','corr(PGA)','P-Value(PGA)','Q-Value(PGA)')],
#            paste0(temp0,'/spearman.xlsx'),rowNames=T)
write.xlsx(corr,paste0(temp0,"/Corr.xlsx"),rowNames = T)
write.xlsx(pvalue,paste0(temp0,"/Pvalue.xlsx"),rowNames = T)
# write.xlsx(corrall,"Corrall.xlsx",rowNames = T)
# write.xlsx(pvalueall,"Pvalueall.xlsx",rowNames = T)


# ###逐步寻优法尝试
# data<-lipiddata[,lipidfilt]
# data$result<-clinical[,'PASI']
# data<-data[!is.na(data$result),]
# mouse.regression<-lm(result~., data=data)
# lm.step<-step(mouse.regression,direction = 'both')
# 
# a=rownames(summary(lm.step)$coefficients)[-1]