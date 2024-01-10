library(openxlsx)
library(ggplot2)
library(readxl)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(cowplot)
library(pROC)
library(dplyr)
library(pals)

setwd("C:/Users/小马很酷/Desktop/To李溢溪/机器学习/10折交叉验证/RA VS HC--随机数种子123随机树500")
data1<-read.xlsx("biomarker.xlsx",sheet = 1,check.names = F,rowNames = T,sep.names = " ")
setwd("C:/Users/小马很酷/Desktop/To李溢溪/机器学习")
raw<-read.xlsx('RA VS HC_RFA.xlsx',rowNames = T,sheet = 2,check.names = T,sep.names = " ")
data3<-as.data.frame(t(raw))#[, -c(ncol(raw),ncol(raw)-1)]
data3$Group<- unlist(sapply(rownames(data3),function(i){strsplit(i,'_')[[1]][1]}))
data3$Group<- unlist(sapply(data3$Group,function(i){ifelse(i=='C','Control',i)}))
data3$Group<-factor(data3$Group,levels = unique(data3$Group))
setwd("C:/Users/小马很酷/Desktop/To李溢溪/机器学习/1000次随机森林/RA VS HC")
data2<-read.xlsx("Gini.xlsx",sheet = 1,check.names = F,sep.names = " ",rowNames = T)

##选前100
data2$average<-sapply(1:nrow(data2),function(i){mean(unlist(data2[i,]))})
data2<-data2[order(data2$average,decreasing = T),]
dataa<-data2[1:100,]
data2<-dataa[order(dataa$average),1:(ncol(dataa)-1)]

data2<-as.data.frame(t(data2))
data3<-data3[,c(rownames(dataa),"Group")]


data<-data2[,1]


for(i in 2:ncol(data2)){
  add<-data2[,i]
  data<-c(data,add)
}

data<-as.data.frame(data)
colnames(data)<-"expr"
data$gene<-rep(colnames(data2),each=1000)
data$Group<-"not hit"
for(j in 1:ncol(data1)){
data[,"Group"][data[,"gene"]==colnames(data1)[j]]="hit"
}

data$gene<-as.factor(data$gene)


# p<-ggboxplot(data,x=reorder("gene","expr"),y="expr",color = "Group",fill = "Group")+##,add = "jitter"
#   scale_color_manual(values = c(hit="red","not hit"="grey"))+
#   scale_fill_manual(values = c(hit="red","not hit"="grey"))+
#   # ggtitle("Accuracy")+
#   theme_bw()+
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     # axis.text = element_text(size = 15,colour = "black"),
#     axis.title = element_text(size = 5,colour = "black"),
#     panel.grid.minor = element_blank(),
#     # legend.position = "none",
#     panel.grid.major = element_blank(),
#     panel.border = element_rect(size = 1.5,colour = "black"),
#     axis.text.y = element_text(size = 5,colour = "black"),
#     axis.text.x = element_text(size = 5,colour = "black",angle = 270),
#   )+
#   coord_flip()+
#   #   compare_means(len ~ supp, data = ToothGrowth, 
#   #                 group.by = "dose", paired = TRUE)
#   # # Box plot facetted by "dose"
#   #   ggpaired(ToothGrowth, x = "supp", y = "len",
#   #               color = "supp", palette = "jama", 
#   #               line.color = "gray", line.size = 0.4,
#   #               facet.by = "dose", short.panel.labs = FALSE)+
#   # stat_compare_means(aes(group = Group),method = "t.test",paired = F)+#, label = "p.signif"
#   labs(x="Name",y="MeanDecreaseAccuracy")

p<-ggplot(data,aes(gene,expr,color = Group,fill = Group))+
  geom_boxplot(aes(x=reorder(gene,expr),y=expr),outlier.size = 0.3,color="black")+
  scale_color_manual(values = c(hit="red","not hit"="grey"))+
  scale_fill_manual(values = c(hit="red","not hit"="grey"))+
  # ggtitle("Accuracy")+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    # axis.text = element_text(size = 15,colour = "black"),
    axis.title = element_text(size = 20,colour = "black"),
    panel.grid.minor = element_blank(),
    # legend.position = "none",
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 1,colour = "black"),
    axis.text.y = element_text(size = 15,colour = "black"),
    axis.text.x = element_text(size = 15,colour = "black",angle = 270),
  )+
  coord_flip()+
  #   compare_means(len ~ supp, data = ToothGrowth, 
  #                 group.by = "dose", paired = TRUE)
  # # Box plot facetted by "dose"
  #   ggpaired(ToothGrowth, x = "supp", y = "len",
  #               color = "supp", palette = "jama", 
  #               line.color = "gray", line.size = 0.4,
  #               facet.by = "dose", short.panel.labs = FALSE)+
  # stat_compare_means(aes(group = Group),method = "t.test",paired = F)+#, label = "p.signif"
  labs(x="Name",y="MeanDecreaseGini")

ggsave("前100MeanDecreaseGini.png",p,width=length(unique(data$gene))*0.12,height = length(unique(data$gene))*0.2,dpi = 600,limitsize = F)

###单独ROC
if (!dir.exists('ROC')){dir.create('ROC')}
auc<-as.data.frame(matrix(nrow = ncol(data3[,-ncol(data3)]),ncol = 3))
colnames(auc)<-c("Genename","aveGini","AUC")
for(i in 1:ncol(data3[,-ncol(data3)])){
  png(paste0("ROC/",colnames(data3)[i],".png"),width = 7,height = 7,units = "in",res = 600)
  proc<-roc(data3$Group,data3[,i],ci=T,main = colnames(data3)[i],plot=T,grid=F,print.auc=T,col="tomato",legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
  pdf(paste0("ROC/",colnames(data3)[i],".pdf"),width = 7,height = 7)
  proc<-roc(data3$Group,data3[,i],ci=T,main = colnames(data3)[i],plot=T,grid=F,print.auc=T,col="tomato",legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
  auc[i,1]<-colnames(data3)[i]
  auc[i,2]<-dataa[colnames(data3)[i],"average"]
  auc[i,3]<-proc$auc
}
write.xlsx(auc,"Gini前100及roc汇总.xlsx")

       