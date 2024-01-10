setwd("D:/马驰宇/桌面/何院专辑/研究院项目/戴主任项目重分析/新数据补值/又有一份新数据")
raw<-read.xlsx("最终统计数据-SLE active-Inactive.xlsx",sheet=1,rowNames=F,check.names=F,sep.names=" ")
data<-raw[,grep("Control|SLE",colnames(raw))]
Control<-data[,grep("Control",data[3,])]
inactive<-data[,grep("Inactive",data[3,])]
active<-data[,grep("Mild|Moderate|Severe",data[3,])]
ksp.Control<-NULL
ksp.inactive<-NULL
ksp.active<-NULL
ksd.Control<-NULL
ksd.inactive<-NULL
ksd.active<-NULL
meanControl<-NULL
meanactive<-NULL
meaninactive<-NULL

for(i in 4:95){
  ksp.Control<-c(ksp.Control,ks.test(as.numeric(unlist(Control[i,]))+runif(length(as.numeric(unlist(Control[i,]))),-0.000005,0.000005),"pnorm",mean(as.numeric(unlist(Control[i,]))),sd(as.numeric(unlist(Control[i,]))))$p.value)
  ksp.active<-c(ksp.active,ks.test(as.numeric(unlist(active[i,]))+runif(length(as.numeric(unlist(active[i,]))),-0.000005,0.000005),"pnorm",mean(as.numeric(unlist(active[i,]))),sd(as.numeric(unlist(active[i,]))))$p.value)
  ksp.inactive<-c(ksp.inactive,ks.test(as.numeric(unlist(inactive[i,]))+runif(length(as.numeric(unlist(inactive[i,]))),-0.000005,0.000005),"pnorm",mean(as.numeric(unlist(inactive[i,]))),sd(as.numeric(unlist(inactive[i,]))))$p.value)
  ksd.Control<-c(ksd.Control,ks.test(as.numeric(unlist(Control[i,]))+runif(length(as.numeric(unlist(Control[i,]))),-0.000005,0.000005),"pnorm",mean(as.numeric(unlist(Control[i,]))),sd(as.numeric(unlist(Control[i,]))))$statistic)
  ksd.active<-c(ksd.active,ks.test(as.numeric(unlist(active[i,]))+runif(length(as.numeric(unlist(active[i,]))),-0.000005,0.000005),"pnorm",mean(as.numeric(unlist(active[i,]))),sd(as.numeric(unlist(active[i,]))))$statistic)
  ksd.inactive<-c(ksd.inactive,ks.test(as.numeric(unlist(inactive[i,]))+runif(length(as.numeric(unlist(inactive[i,]))),-0.000005,0.000005),"pnorm",mean(as.numeric(unlist(inactive[i,]))),sd(as.numeric(unlist(inactive[i,]))))$statistic)
  meanControl<-c(meanControl,mean(as.numeric(unlist(Control[i,]))))
  meanactive<-c(meanactive,mean(as.numeric(unlist(active[i,]))))
  meaninactive<-c(meaninactive,mean(as.numeric(unlist(inactive[i,]))))
}
pvalueca<-NULL
pvalueci<-NULL
pvalueia<-NULL
qvalueca<-NULL
qvalueci<-NULL
qvalueia<-NULL
for(i in 4:95){
  if(ksp.inactive[i-3]>0.05&ksp.active[i-3]>0.05){
    varia<-var.test(as.numeric(unlist(active[i,],use.names = F)),as.numeric(unlist(inactive[i,],use.names = F)))$p.value > 0.05
    pvalueia<-c(pvalueia,t.test(as.numeric(unlist(active[i,],use.names = F)), as.numeric(unlist(inactive[i,],use.names = F)), var.equal = varia)$p.value)
  }else{
    pvalueia<-c(pvalueia,wilcox.test(as.numeric(unlist(active[i,])),as.numeric(unlist(inactive[i,])))$p.value)
  }
  if(ksp.active[i-3]>0.05&ksp.Control[i-3]>0.05){
    varca<-var.test(as.numeric(unlist(Control[i,],use.names = F)),as.numeric(unlist(active[i,],use.names = F)))$p.value > 0.05
    pvalueca<-c(pvalueca,t.test(as.numeric(unlist(Control[i,],use.names = F)), as.numeric(unlist(active[i,],use.names = F)), var.equal = varca)$p.value)
  }else{
    pvalueca<-c(pvalueca,wilcox.test(as.numeric(unlist(Control[i,])),as.numeric(unlist(active[i,])))$p.value)
  }
  if(ksp.inactive[i-3]>0.05&ksp.Control[i-3]>0.05){
    varci<-var.test(as.numeric(unlist(Control[i,],use.names = F)),as.numeric(unlist(inactive[i,],use.names = F)))$p.value > 0.05
    pvalueci<-c(pvalueci,t.test(as.numeric(unlist(Control[i,],use.names = F)), as.numeric(unlist(inactive[i,],use.names = F)), var.equal = varci)$p.value)
  }else{
    pvalueci<-c(pvalueci,wilcox.test(as.numeric(unlist(Control[i,])),as.numeric(unlist(inactive[i,])))$p.value)
  }
}
qvalueca<-p.adjust(pvalueca,method = "BH")
qvalueci<-p.adjust(pvalueci,method = "BH")
qvalueia<-p.adjust(pvalueia,method = "BH")
output<-cbind(raw[,1:3],data)
output$meanControl<-c(NA,NA,NA,meanControl)
output$meanactive<-c(NA,NA,NA,meanactive)
output$meaninactive<-c(NA,NA,NA,meaninactive)
output$ks.pvalue.Control<-c(NA,NA,NA,ksp.Control)
output$ks.statistic.Control<-c(NA,NA,NA,ksd.Control)
output$ks.pvalue.inactive<-c(NA,NA,NA,ksp.inactive)
output$ks.statistic.inactive<-c(NA,NA,NA,ksd.inactive)
output$ks.pvalue.active<-c(NA,NA,NA,ksp.active)
output$ks.statistic.active<-c(NA,NA,NA,ksd.active)
output$fc.active.Control<-c(NA,NA,NA,meanactive/meanControl)
output$fc.inactive.Control<-c(NA,NA,NA,meaninactive/meanControl)
output$fc.active.inactive<-c(NA,NA,NA,meanactive/meaninactive)
output$pvalue.active.inactive<-c(NA,NA,NA,pvalueia)
output$qvalue.active.inactive<-c(NA,NA,NA,qvalueia)
output$pvalue.Control.inactive<-c(NA,NA,NA,pvalueci)
output$qvalue.Control.inactive<-c(NA,NA,NA,qvalueci)
output$pvalue.active.Control<-c(NA,NA,NA,pvalueca)
output$qvalue.active.Control<-c(NA,NA,NA,qvalueca)
write.xlsx(output,"显著性检验.xlsx",rowNames=F)

library(openxlsx)
library(ggplot2)
library(readxl)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(reshape2)
setwd("D:/马驰宇/桌面/何院专辑/研究院项目/戴主任项目重分析/新数据补值/又有一份新数据")
raw<-read.xlsx("显著性检验Control-active-inactive.xlsx",sheet = 1,rowNames=F,sep.names=" ")
for(i in 110:230){
  colnames(raw)[i]<-str_replace_all(colnames(raw)[i],"SLE",raw[3,i])
}
raw<-raw[-c(1:3),]
rownames(raw)<-raw[,1]
raw<-raw[,-c(1:3)]
for(i in colnames(raw)){
  raw[,i]<-as.numeric(raw[,i])
}

temps<-'Control-active'
if(!dir.exists(temps)){dir.create(temps)}
test<-NULL
for(i in 1:nrow(raw)){
  if(as.numeric(raw$ks.pvalue.Control[i])>0.05&as.numeric(raw$ks.pvalue.active[i])>0.05){
    test<-c(test,"t.test")
  }else{
    test<-c(test,"wilcox.test")
  }
}
data<-raw[,grep("Control|Mild|Severe|Moderate",colnames(raw)[1:227])]
data<-as.data.frame(t(data))
data$group<-str_remove_all(rownames(data),"\\d")
data$group[data$group=="Mild"|data$group=="Severe"|data$group=="Moderate"]<-"Active"
for(i in 1:92){
  data[,i]<-as.numeric(data[,i])
  p<-ggplot(data,aes(group,data[,i]),size=0.4)+
    stat_boxplot(geom = "errorbar",size=0.4,width = 0.2,aes(color = group))+
    geom_boxplot(aes(color = group),size = 0.4,fill = "white",width = 0.4,outlier.alpha = 0.8)+
    # geom_jitter(aes(color = group),width =0.175,shape = 16,size=2)+
    # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,aes(color = group,fill = group))+
    scale_color_manual(values = c(Control="#9B30FF",Active="#FF8247"))+
    # scale_fill_manual(values = c(Control="#9B30FF",SLE="#FF8247"))+
    ggtitle(colnames(data)[i])+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 15,colour = "black"),
      axis.title = element_text(size = 15,colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.border = element_rect(size = 1.5,colour = "black")
    )+
    #   compare_means(len ~ supp, data = ToothGrowth, 
    #                 group.by = "dose", paired = TRUE)
    # # data plot facetted by "dose"
    #   ggpaired(ToothGrowth, x = "supp", y = "len",
    #               color = "supp", palette = "jama", 
    #               line.color = "gray", line.size = 0.4,
    #               facet.by = "dose", short.panel.labs = FALSE)+
    stat_compare_means(comparisons = list(c(data$group[1],data$group[nrow(data)])), method = test[i], paired = F,label = "p.signif")+
    # stat_compare_means()+
    labs(x=NULL,y="Expression")+
    coord_cartesian(ylim = c(0,1.1*max(data[,i])))+
    scale_y_continuous(expand = c(0, 0))
  ggsave(paste0(temps,"/",colnames(data)[i],"-",test[i],"-",raw$pvalue.active.Control[i],".png"), p, width = 3.2, height = 4.8,units = "in",dpi = 600)
  ggsave(paste0(temps,"/",colnames(data)[i],"-",test[i],"-",raw$pvalue.active.Control[i],".pdf"), p, width = 3.2, height = 4.8)
}