library(openxlsx)
library(ggplot2)
library(readxl)
library(ggpubr)
library(ggsignif)
library(ggsci)
library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(reshape2)

###一个一个画
for(i in colnames(datas)[-11]){
  x<-datas$group
  y<-datas[,i]
  yname<-paste0(rownames(datas),i)
  yname<-str_replace_all(yname,paste0("1",i),i)
  data<-data.frame(x=x,y=y,yname=yname)
  
  p<-ggplot(data,aes(x,y),size=0.4)+
    stat_boxplot(geom = "errorbar",size=0.4,width = 0.2,aes(color = x))+
    geom_boxplot(aes(color = x),position = position_dodge(width = 0.2),notch=TRUE,notchwidth=0.8,size = 0.4,fill = "white",width = 0.4,outlier.alpha = 0)+
    # geom_jitter(aes(color = Group),width =0.175,shape = 16,size=2)+
    geom_point(size=1,aes(color = x,fill = x))+
    geom_hline(yintercept = 0,color = "grey",linetype = "dashed")+
    geom_line(aes(group=yname) ,size=0.5,color="grey50")+
    scale_color_manual(values = c(HC="#9B30FF",RT="#FF8247"))+
    # scale_fill_manual(values = c(Control="#9B30FF",SLE="#FF8247"))+
    ggtitle(i)+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 10,colour = "black"),
      axis.title = element_text(size = 10,colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.border = element_rect(size = 1,colour = "black")
    )+
    #   compare_means(len ~ supp, data = ToothGrowth, 
    #                 group.by = "dose", paired = TRUE)
    # # Box plot facetted by "dose"
    #   ggpaired(ToothGrowth, x = "supp", y = "len",
    #               color = "supp", palette = "jama", 
    #               line.color = "gray", line.size = 0.4,
    #               facet.by = "dose", short.panel.labs = FALSE)+
    stat_compare_means(comparisons = list(c(data$x[1],data$x[nrow(datas)])), method = "t.test", paired = F,label = "p.signif")+
    # stat_compare_means()+
    labs(x=NULL,y="Correlation coefficient")
  # coord_cartesian(ylim = c(0,1.1*max(datas[,i])))
  # scale_y_continuous(expand = c(0, 0))
  ggsave(paste0(i,".png"), p, width = 3.2, height = 4.8,units = "in",dpi = 600)
  ggsave(paste0(i,".pdf"), p, width = 3.2, height = 4.8)
}




###全部一起

for(sheet in excel_sheets("SLE_PBMCVSNC_PBMC_Gene_differential_expression(已自动还原).xlsx")){
  if(!dir.exists(sheet)){dir.create(sheet)}
  raw<-read.xlsx("SLE_PBMCVSNC_PBMC_Gene_differential_expression(已自动还原).xlsx",sheet = sheet,colNames = T,rowNames = T,check.names = F,sep.names = " ")
  rownames(raw)<-raw[,1]
  raw<-raw[,-1]
  data<-raw[,c(7:12)]
  # data<-raw
  data<-as.data.frame(t(data))
  expr<-NULL
  for(i in 1:ncol(data)){
    expr<-c(expr,data[,i])}
  gene<-rep(colnames(data),each = nrow(data))
  Group<-rep(rownames(data),times = ncol(data))
  box<-data.frame(expr,gene,Group)
  box$Group<-sapply(sapply(rownames(data),function(ii){strsplit(ii,"_")[[1]][1]}),function(iii){strsplit(iii,"FPKM.")[[1]][2]})
  box$Group<-ifelse(box$Group=="NC","Control","SLE")
  # box$Group<-sapply(rownames(data),function(ii){strsplit(ii,"_")[[1]][1]})
    # box$Group<-factor(box$Group,levels = box$Group)
  p<-ggboxplot(data,x="x",y="y",color = "color",add = "jitter")+##,add = "jitter"
    scale_color_manual(values = c(up="#FF3030",down="#436EEE"))+
    # ggtitle(sheet)+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      # axis.text = element_text(size = 15,colour = "black"),
      axis.title = element_text(size = 15,colour = "black"),
      panel.grid.minor = element_blank(),
      # legend.position = "none",
      panel.grid.major = element_blank(),
      panel.border = element_rect(size = 1.5,colour = "black"),
      axis.text.y = element_text(size = 15,colour = "black"),
      axis.text.x = element_text(size = 15,colour = "black",angle = 0),
    )+
      #   compare_means(len ~ supp, data = ToothGrowth, 
      #                 group.by = "dose", paired = TRUE)
      # # Box plot facetted by "dose"
      #   ggpaired(ToothGrowth, x = "supp", y = "len",
      #               color = "supp", palette = "jama", 
      #               line.color = "gray", line.size = 0.4,
      #               facet.by = "dose", short.panel.labs = FALSE)+
    # stat_compare_means(aes(group = color),method = "t.test",paired = F)+#, label = "p.signif"
    labs(x=NULL,y="Correlation")
  ggsave(paste0(sheet,"/",colnames(box)[1],".png"), p, width = 9.6, height = 4.8,units = "in",dpi = 600)
  ggsave(paste0(sheet,"/all-1.pdf"), p, width = 9.6, height = 4.8)
  write.xlsx(box,paste0(sheet,"/all.xlsx"))
}


####散点
data<-read.xlsx("总结.xlsx",sheet = 3)
data$y<-abs(data$y)
p <- ggplot(data, aes(x=x, y=y)) + 
  geom_jitter(aes(colour=Correlation),shape=18,width = 0.4,size=4) + ##geom_quasirandom
  scale_color_manual(values = c(Positive="#FF3030",Negative="#436EEE"))+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.key=element_blank(),
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 15),
                     axis.ticks = element_blank(),
                     legend.title = element_text(size = 15),
                     legend.text = element_text(size = 15),
                     legend.position = c(0.85,0.8))+
  coord_cartesian(ylim = c(0.1,0.55))+
  labs(x = NULL, y = "Correlation")
  # theme(legend.position="none")
ggsave("jitterplot2.png",p, width = 7.2, height = 4.8,units = "in",dpi = 600)
ggsave("jitterplot2.pdf",p, width = 7.2, height = 4.8)


raw<-read.xlsx("data.xlsx",sheet = 4,rowNames = T)
raw$x<-factor(raw$x,levels = unique(raw$x))
color<-jet(length(unique(raw$x)))
names(color)<-unique(raw$x)

###箱线图叠加折线图
for(i in colnames(raw)[-10]){
  data<-data.frame(x=raw$x,y=raw[,i])
  data %>% 
    group_by(x) %>% 
    summarise(average=mean(y))->dat
  p<-ggplot(data,aes(x,y),size=0.4)+
    stat_boxplot(geom = "errorbar",size=0.4,width = 0.2,aes(color = x))+
    geom_boxplot(aes(color = x),position = position_dodge(width = 0.2),size = 0.4,fill = "white",width = 0.4,outlier.alpha = 0)+
    # geom_jitter(aes(color = Group),width =0.175,shape = 16,size=2)+
    geom_point(size=1,aes(color = x))+
    # geom_hline(yintercept = 0,color = "grey",linetype = "dashed")+
    # geom_line(aes(group=x) ,size=0.5,color="grey50")+
    scale_color_manual(values = color)+
    geom_line(data=dat,aes(x=x,y=average,group=1),color="black")+
    # scale_fill_manual(values = c(Control="#9B30FF",SLE="#FF8247"))+
    ggtitle(i)+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 10,colour = "black"),
      axis.title = element_text(size = 10,colour = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.border = element_rect(size = 1,colour = "black")
    )+
    #   compare_means(len ~ supp, data = ToothGrowth, 
    #                 group.by = "dose", paired = TRUE)
    # # Box plot facetted by "dose"
    #   ggpaired(ToothGrowth, x = "supp", y = "len",
    #               color = "supp", palette = "jama", 
    #               line.color = "gray", line.size = 0.4,
    #               facet.by = "dose", short.panel.labs = FALSE)+
    stat_compare_means(comparisons = list(c("db_m","db_db"),c("db_db","Empa")), method = "t.test", paired = F,label = "p.signif")+
    # stat_compare_means()+
    labs(x=NULL,y=NULL)
  # coord_cartesian(ylim = c(0,1.1*max(datas[,i])))
  # scale_y_continuous(expand = c(0, 0))
  ggsave(paste0(i,".png"), p, width = 6, height = 4.8,units = "in",dpi = 600)
  ggsave(paste0(i,".pdf"), p, width = 6, height = 4.8)
}
