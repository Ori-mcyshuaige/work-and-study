library(ropls)
library(ggplot2)
library(data.table)###多线程读文件
library(pals)
library(stringr)
library(psych)####相关性分析
library(openxlsx)
options(stringsAsFactors = F,timeout = 200)






ModelScorePlot <- function(model, class.id, output,colanno,ci=0.95,ellipse=F,labelshow=F) {
  # PCA/OPLS-DA模型得分图
  # 
  if(model@typeC == 'PCA' | model@typeC == 'PLS-DA') {
    df.p <- as.data.frame(model@scoreMN[, 1:2])
  } else if(model@typeC == 'OPLS-DA') {
    df.p <- cbind(
      as.data.frame(model@scoreMN[, 1]), 
      as.data.frame(model@orthoScoreMN[,1])
    )
  }
  colnames(df.p) <- c('x', 'y')
  # write.xlsx(df.p,file = paste0(output, '/', model@typeC, 'data.xlsx'),row.names=TRUE)
  n <- nrow(df.p)
  hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(ci, 2, (n-2))
  rv <- seq(0, 2*pi, length.out = 100)
  df.ell <- data.frame(  # 置信区间数据
    x = sqrt(var(df.p$x)*hfn)*cos(rv), 
    y = sqrt(var(df.p$y)*hfn)*sin(rv)
  )
  p <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_polygon(data = df.ell, aes(x, y), color = 'black', fill = NA)
  if(ellipse) {
    p <- p + stat_ellipse(
      data = df.p, geom = 'polygon', level = ci, 
      aes(x, y, fill = class.id, color = class.id), alpha = I(0.1)
    ) +
      stat_ellipse(
        data = df.p, geom = 'blank', level = ci, 
        aes(-x, -y, fill = class.id)
      )
  }
  if(labelshow) {
    p <- p + geom_point(
      data = df.p, size = 6, aes(x, y, shape = class.id, color = class.id)
    ) +
      geom_blank(
        data = df.p, aes(-x, -y, shape = class.id, color = class.id)
      ) +
      scale_shape_manual(  # 形状参数
        values = rep_len(
          #c("\u2605","\u25C4","\u25BC","\u25B2"), 
          c(16, 15, 17, 18),
          length(levels(class.id))
        )[sort(unique(as.numeric(class.id)))]
      ) +
      scale_color_manual(  # 颜色参数
        values = colanno
      ) +
      scale_fill_manual(
        values = colanno
      ) +
      geom_label_repel(
        data = df.p, 
        mapping = aes(x, y, label = rownames(df.p)), 
        color = 'black', size = 3, 
        label.padding = unit(0.2, 'lines'), 
        point.padding = unit(0.5, 'lines'), 
        min.segment.length = unit(0.1, "lines"), 
        segment.color = 'grey50', segment.size = 1, 
        show.legend = F
      ) +
      theme_bw() +
      theme(
        legend.title = element_blank(), 
        legend.key = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 1)
        # axis.text = element_text(size = 15, color = "black",
        #                          face = "plain", vjust = 0.5, hjust = 0.5),
        # axis.title = element_text(size = 15, color = "black", 
        #                           face = "plain", vjust = 0.5, hjust = 0.5)
      ) +
      
      if(model@typeC == 'PCA') {
        labs(x = 'PC[1]', y = 'PC[2]')
      } else if(model@typeC == 'OPLS-DA') {
        labs(x = 't[1]P', y = 't[1]O')
      } else if(model@typeC == 'PLS-DA') {
        labs(x = 't[1]', y = 't[2]')
      }
  }
  if(labelshow == F) {
    p <- p + geom_point(
      data = df.p, size = 5, aes(x, y, shape = class.id, color = class.id)
    ) +
      geom_blank(
        data = df.p, aes(-x, -y, shape = class.id, color = class.id)
      ) +
      scale_shape_manual(  # 形状参数
        values = rep_len(
          #c("\u2605","\u25C4","\u25BC","\u25B2"), 
          c(16, 15, 17, 18),
          length(levels(class.id))
        )[sort(unique(as.numeric(class.id)))]
      ) +
      scale_color_manual(  # 颜色参数
        values = colanno
      ) +
      scale_fill_manual(
        values = colanno
      ) +
      theme_bw() +
      theme(
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 1)
        # axis.text = element_text(size = 15, color = "black",
        #                          face = "plain", vjust = 0.5, hjust = 0.5),
        # axis.title = element_text(size = 15, color = "black", 
        #                           face = "plain", vjust = 0.5, hjust = 0.5)
      ) +
      if(model@typeC == 'PCA') {
        labs(x = paste0('PC1(',round(100*model@pcaVarVn['p1']/sum(model@pcaVarVn),2),'%)'),
             y = paste0('PC2(',round(100*model@pcaVarVn['p2']/sum(model@pcaVarVn),2),'%)'))
      } else if(model@typeC == 'OPLS-DA') {
        labs(x = 't[1]P', y = 't[1]O')
      } else if(model@typeC == 'PLS-DA') {
        labs(x = 't[1]', y = 't[2]')
      }
  }
  ggsave(  # 输出位置
    paste0(output, '/', model@typeC, ' score plot.jpg'), p, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', model@typeC, ' score plot.pdf'), p,
    device = cairo_pdf,
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
}

ModelScorePlot3DNew <- function(model, class.id, output,colanno) {
  # 新的3D PCA得分图~~撒花~
  # class.id参数必须是因子，不然会出现图例的排序问题
  # 
  require(plotly)
  require(htmlwidgets)
  df.data <- as.data.frame(model@scoreMN[, 1:3])
  df.data <- cbind(class.id, df.data)
  colnames(df.data) <- c("class.id", "PC1", "PC2", "PC3")
  # write.xlsx(df.data,file = paste0(output, '/', model@typeC, '_3D_data.xlsx'),row.names=TRUE)
  axis.x <- list(
    range = c(- max(abs(df.data$PC1)), max(abs(df.data$PC1)))
  )
  axis.y <- list(
    range = c(- max(abs(df.data$PC2)), max(abs(df.data$PC2)))
  )
  axis.z <- list(
    range = c(- max(abs(df.data$PC3)), max(abs(df.data$PC3)))
  )
  p <- plot_ly(
    df.data, x = ~PC1, y = ~PC2, z = ~PC3, 
    type = "scatter3d", mode = "markers", 
    color = ~class.id, 
    colors = colanno, 
    symbol = ~class.id, 
    symbols = rep_len(
      #c("\u2605","\u25C4","\u25BC","\u25B2"),
      c(16, 15, 17, 18),
      length.out = length(levels(class.id))
    ), 
    text = ~row.names(df.data), 
    marker = list(size = 8, opacity = 0.8)
  ) %>% 
    layout(
      legend = list(
        bgcolor = "#FFFFFF",
        bordercolor = "#000000"
        
      ), 
      scene = list(
        xaxis = axis.x, yaxis = axis.y, zaxis = axis.z
      )
    )
  old.wd <- getwd()
  setwd(output)
  on.exit(setwd(old.wd))
  htmlwidgets::saveWidget(p, file = paste0(model@typeC," score plot 3D.html"))
}

Splot<-function(df.model,model,output,scaleC){
  plotdata<-log(df.model[,colnames(df.model) %in% names(model@vipVn)],base = 10)

  library(parallel)
  x <- 1:ncol(plotdata)
  cl <- makeCluster(8) # 初始化8核心集群
 
  results <- parLapply(cl,x,function(i,plotdata,scaleC){
    if(scaleC=='center'){
      s<-data.frame(plotdata[,i]-mean(plotdata[,i]))
      colnames(s)<-colnames(plotdata)[i]
      s
    }else if(scaleC=='standard'){
      s<-data.frame((plotdata[,i]-mean(plotdata[,i]))/sd(plotdata[,i]))
      colnames(s)<-colnames(plotdata)[i]
      s
    }else{
      s<-data.frame((plotdata[,i]-mean(plotdata[,i]))/sqrt(sd(plotdata[,i])))
      colnames(s)<-colnames(plotdata)[i]
      s
    }
  },plotdata,scaleC) # lapply的并行版本
  jinghua <- data.frame(results) # 整合结果
  stopCluster(cl) # 关闭集群
  
  x=cov(model@scoreMN[,1],jinghua)
  y=corr.test(model@scoreMN[,1],jinghua)
  df.p<-data.frame(x=t(x),y=t(y$r),row.names = colnames(jinghua))
  write.xlsx(df.p,file = paste0(output, '/', model@typeC, '_splot.xlsx'),row.names=TRUE)
  
  p <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(
      data = df.p, size = 1, aes(x, y)
    ) +
    geom_blank(
      data = df.p, aes(-x, -y)
    )+
    labs(x = 'Cov(tp,X)', y = 'Corr(tp,X)') +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 1)
    ) 
      
      
  
  ggsave(  # 输出位置
    paste0(output, '/', model@typeC, ' S-plot.jpg'), p, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', model@typeC, ' S-plot.pdf'), p,
    device = cairo_pdf,
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  
}



Vplot<-function(df.model,model,output,scaleC){
  plotdata<-log(df.model[,colnames(df.model) %in% names(model@vipVn)],base = 10)
  
  library(parallel)
  x <- 1:ncol(plotdata)
  cl <- makeCluster(8) # 初始化8核心集群
  
  results <- parLapply(cl,x,function(i,plotdata,scaleC){
    if(scaleC=='center'){
      s<-data.frame(plotdata[,i]-mean(plotdata[,i]))
      colnames(s)<-colnames(plotdata)[i]
      s
    }else if(scaleC=='standard'){
      s<-data.frame((plotdata[,i]-mean(plotdata[,i]))/sd(plotdata[,i]))
      colnames(s)<-colnames(plotdata)[i]
      s
    }else{
      s<-data.frame((plotdata[,i]-mean(plotdata[,i]))/sqrt(sd(plotdata[,i])))
      colnames(s)<-colnames(plotdata)[i]
      s
    }
  },plotdata,scaleC) # lapply的并行版本
  jinghua <- data.frame(results) # 整合结果
  stopCluster(cl) # 关闭集群
  
  # x=cov(model@scoreMN[,1],jinghua)
  x=corr.test(model@scoreMN[,1],jinghua)
  df.p<-data.frame(y=(model@vipVn),x=t(x$r),row.names = colnames(jinghua))
  write.xlsx(df.p,file = paste0(output, '/', model@typeC, '_Vplot.xlsx'),row.names=TRUE)
  
  p <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(
      data = df.p, size = 1, aes(x, y)
    ) +
    geom_blank(
      data = df.p, aes(-x, 0)
    )+
    labs(y = 'VIP', x = 'Corr(tp,X)') +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = 'black', size = 1)
    ) 
  
  
  
  ggsave(  # 输出位置
    paste0(output, '/', model@typeC, ' V-plot.jpg'), p, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', model@typeC, ' V-plot.pdf'), p,
    device = cairo_pdf,
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  
}





getpca<-function(data,class.id,output,colanno=NA,log10L=T,scaleC='center',ci=0.95,ellipse=F,labelshow=F){
  df.model <- data
  pca <- opls(
    x = df.model, 
    log10L = log10L, scaleC = scaleC, 
    crossvalI = min(nrow(df.model), 7), 
    printL = F, plotL = F
  )
  if(pca@summaryDF$pre < 3) {  # 主成分不少于3，方便做三维图
    pca <- opls(
      x = df.model, predI = 3, 
      log10L = log10L, scaleC = scaleC, 
      crossvalI = min(nrow(df.model), 7), 
      printL = F, plotL = F
    )
  }
  temps<-str_split(output,'/')[[1]]
  for (ii in 1:length(temps)) {
    temp0<-paste0(temps[1:ii],collapse = '/')
    if(!dir.exists(temp0)){dir.create(temp0)}
  }
  
  if(is.na(colanno)){
    colanno<-jet(length(unique(class.id)))
    names(colanno)<-unique(class.id)
  }
  ModelScorePlot(  # 2D得分图
    model = pca, class.id = class.id, colanno=colanno,
    output = output,ci=ci,ellipse=ellipse,labelshow=labelshow
  )
  ModelScorePlot3DNew(  # 3D得分图
    model = pca, class.id = class.id, colanno=colanno,
    output = output
  )
  
  pcdata<-pca@scoreMN
  write.xlsx(pcdata,file = paste0(output, '/', pca@typeC, 'data.xlsx'),row.names=TRUE)
  write.csv(pca@pcaVarVn,file = paste0(output, '/', pca@typeC, '_value.csv'),row.names=TRUE)
  write.xlsx(pca@loadingMN,file = paste0(output, '/', pca@typeC, '_loadingMN.xlsx'),row.names=TRUE)
  statdata<- cbind(pca@xMeanVn,pca@xSdVn)
  colnames(statdata)<-c('xMeanVn','xSdVn')
  write.xlsx(statdata,file = paste0(output, '/', pca@typeC, '_xVn.xlsx'),row.names=TRUE)
}


getoplsdata<-function(data,class.id,output,colanno=NA,log10L=T,scaleC='center',ci=0.95,ellipse=F,labelshow=F){
  df.model <- data
  set.seed(123)
  opls.da <- opls(
    x = df.model, y=as.character(class.id),
    predI = 1, ortho = 1, permI = 200, 
    log10L = log10L, scaleC = scaleC, 
    crossvalI = min(nrow(df.model), 7), 
    printL = F, plotL = F
  )
 
  temps<-str_split(output,'/')[[1]]
  for (ii in 1:length(temps)) {
    temp0<-paste0(temps[1:ii],collapse = '/')
    if(!dir.exists(temp0)){dir.create(temp0)}
  }
  
  if(is.na(colanno)){
    colanno<-jet(length(unique(class.id)))
    names(colanno)<-unique(class.id)
  }
  ModelScorePlot(  # 2D得分图
    model = opls.da, class.id = class.id, colanno=colanno,
    output = output,ci=ci,ellipse=ellipse,labelshow=labelshow
  )
  Splot(df.model,opls.da,output = output,scaleC = scaleC)
  Vplot(df.model,opls.da,output = output,scaleC = scaleC)
  # if(scaleC!='standard'){
  #   Splot(df.model,opls.da,output = output,scaleC = scaleC)
  # }else{
  #   Vplot(df.model,opls.da,output = output,scaleC = scaleC)
  # }
  
 
  

 
  statdata<- cbind(opls.da@vipVn,opls.da@orthoVipVn,opls.da@coefficientMN,opls.da@loadingMN,opls.da@weightMN,opls.da@orthoLoadingMN,opls.da@orthoWeightMN)
  colnames(statdata)<-c('vipVn','orthoVipVn','coefficientMN','loadingMN','weightMN','orthoLoadingMN','orthoWeightMN')
  write.xlsx(statdata,file = paste0(output, '/', opls.da@typeC, '_value.xlsx'),row.names=TRUE)
  statdata<- cbind(opls.da@xMeanVn,opls.da@xSdVn)
  colnames(statdata)<-c('xMeanVn','xSdVn')
  write.xlsx(statdata,file = paste0(output, '/', opls.da@typeC, '_xVn.xlsx'),row.names=TRUE)
  statdata<- cbind(opls.da@scoreMN,opls.da@orthoScoreMN)
  colnames(statdata)<-c('scoreMN','orthoScoreMN')
  write.xlsx(statdata,file = paste0(output, '/', opls.da@typeC, '.xlsx'),row.names=TRUE)
}








posraw<-fread('../LCQEPOSrsd filtered data.csv',stringsAsFactors = F,nThread=8)
posIS<-which(posraw$`MS2 name`=='IS')
posnor<-posraw
posnor<-data.frame(posnor)
posmin<-min(posnor[,7:ncol(posnor)],na.rm = T)
posnor[,7:ncol(posnor)][is.na(posnor[,7:ncol(posnor)])]<-posmin/2
for (jj in 7:ncol(posnor)) {
  posnor[,jj]<-as.numeric(posnor[,jj])/as.numeric(posnor[posIS,jj])
}
# rownames(posnor)<-unlist(sapply(posnor$id, function(i){paste0(i,'_POS')}))
rownames(posnor)<-unlist(sapply(1:nrow(posnor), function(i){paste0(i,'_POS')}))



craw<-read.table('正常人分组.txt',header=T)
cgroup<-craw$Group
names(cgroup)<-craw$Sample
negraw<-fread('../LCQENEGrsd filtered data.csv',stringsAsFactors = F,nThread=8)
negIS<-which(negraw$`MS2 name`=='IS')
negnor<-negraw
negnor<-data.frame(negnor)
negmin<-min(negnor[,7:ncol(negnor)],na.rm = T)
negnor[,7:ncol(negnor)][is.na(negnor[,7:ncol(negnor)])]<-posmin/2
for (jj in 7:ncol(negnor)) {
  negnor[,jj]<-as.numeric(negnor[,jj])/as.numeric(negnor[negIS,jj])
}
# rownames(negnor)<-unlist(sapply(negnor$id, function(i){paste0(i,'_NEG')}))
rownames(negnor)<-unlist(sapply(1:nrow(negraw), function(i){paste0(i,'_NEG')}))

negnor<-negnor[,colnames(posnor)]



totalnor<-rbind(posnor,negnor)

pcadata<-t(totalnor[,7:ncol(totalnor)])
groupdata<-unlist(sapply(rownames(pcadata), function(i){ifelse(str_starts(i,'QC'),'QC',str_split(i,'\\.')[[1]][1])}))


getpca(pcadata,class.id = factor(groupdata,levels = unique(groupdata)),scaleC = 'standard',
       output = 'PCA')
# getpca(t(posnor[,7:ncol(posnor)]),class.id = factor(groupdata,levels = unique(groupdata)),scaleC = 'standard',
#        output = 'POS')
# getpca(t(negnor[,7:ncol(negnor)]),class.id = factor(groupdata,levels = unique(groupdata)),scaleC = 'standard',
#        output = 'NEG')

getpca(t(posnor[,names(cgroup)]),class.id = factor(cgroup,levels = unique(cgroup)),scaleC = 'standard',
       output = 'POS_C')

getoplsdata(t(posnor[,names(cgroup)]),class.id = factor(cgroup,levels = unique(cgroup)),scaleC = 'standard',
            output = 'POS_C/UV')
getoplsdata(t(posnor[,names(cgroup)]),class.id = factor(cgroup,levels = unique(cgroup)),scaleC = 'pareto',
            output = 'POS_C/Par')
getoplsdata(t(posnor[,names(cgroup)]),class.id = factor(cgroup,levels = unique(cgroup)),#scaleC = 'pareto',
            output = 'POS_C/Ctr')
# set.seed(123)  # 固定随机数种子
# pls.da <- opls(
#   x = t(posnor[,names(cgroup)]), y = cgroup, 
#   permI = 200, log10L = T, 
#   scaleC = 'standard', crossvalI = min(nrow(t(posnor[,names(cgroup)])), 7), 
#   printL = F, plotL = F
# )
# set.seed(123)  # 固定随机数种子
# opls.da <- opls(
#   x = t(posnor[,names(cgroup)]), y = cgroup, 
#   predI = 1, ortho = 1, permI = 200, 
#   log10L = T, scaleC = 'standard', 
#   crossvalI = min(nrow(t(posnor[,names(cgroup)])), 7), 
#   printL = F, plotL = F
# )

