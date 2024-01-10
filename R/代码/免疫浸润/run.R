install.packages('e1071')

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("preprocessCore")


inputFile="normalize.txt"      #?????????ļ?
setwd("C:\\biowolf\\geoCRG\\11.CIBERSORT")      #???ù???Ŀ¼
source("CIBERSORT.R")       #???ð?

#????ϸ??????????
outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)

#?????߽??????????ˣ????ұ???????ϸ??????????
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)

###画图
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
raw <- read.table('CIBERSORT-Results.txt', sep="\t",check.names = F,row.names = 1,header = T)
data <- data.frame(sample=rep(rownames(raw),ncol(raw)),
                   proportion=unlist(sapply(colnames(raw),function(x){rbind(as.data.frame(raw[,x]))})),
                   celltype=rep(colnames(raw),each=nrow(raw)))
data$sample <- sapply(data$sample,function(x){paste0(strsplit(x,'_')[[1]][4],'_',strsplit(x,'_')[[1]][3])})

data$group <- sapply(data$sample, function(x){strsplit(x,'_')[[1]][1]})

# data$sample <- factor(data$sample,levels=unique(data$sample))
# data$group <- factor(data$group,levels=data$group)

getcolorname_scale_color_gradient<-function(x,low='#00007F',high='#7F0000', space = "Lab"){
  require(scales)
  get_palette<-seq_gradient_pal(low,high, space)
  get_rescaler<-function(x=x,to=c(0,1),from=range(x,na.rm=TRUE)){rescale(x,to,from)}
  return(get_palette(get_rescaler(x)))
}
data$color <- getcolorname_scale_color_gradient(as.numeric(factor(data$group))-1)
data <- rbind(data[data$group=='control',],data[data$group!='control',])
write.csv(data,'boxplotdata.csv')
p <- ggplot(data,aes(sample,proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity",just = 0.5,width = 1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "white",size = 1),
    legend.position = 'bottom',
    axis.ticks.y = element_line(size = .1,linetype=1),
    axis.ticks.length.x = unit(0.5,'cm'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(size = 2.1,linetype=2,color = rep(unique(data$color),each=64))
    # axis.text.x= element_text(face = "bold", color = 'black', size = 12,angle=90)
  )+
  # guides(fill=guide_legend(ncol=1))+
  labs(x='',
       y='Relative Percent',
       fill='')+
  scale_y_continuous(limits = c(0,1.001),expand = c(0,0))+
  guides(col = guide_legend(ncol= 1))
  # facet_grid(~data$group,scales= "free",space= "free")
ggsave('./percentbar.pdf',p,width = 9,height = 5,dpi = 600)
ggsave('./percentbar.png',p,width = 9,height = 5,dpi = 600)

p<-ggboxplot(data,x="celltype",y="proportion",color = "group",add = "jitter")+##,add = "jitter"
  scale_color_manual(values = c(MDD="#FF3030",control="#436EEE"))+
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
    axis.text.y = element_text(size = 12,colour = "black"),
    axis.text.x = element_text(size = 12,colour = "black",angle = 45,vjust = 1,hjust = 1),
  )+
    # compare_means(len ~ supp, data = ToothGrowth,
    #               group.by = "dose", paired = TRUE),
  # # Box plot facetted by "dose"
  #   ggpaired(ToothGrowth, x = "supp", y = "len",
  #               color = "supp", palette = "jama", 
  #               line.color = "gray", line.size = 0.4,
  #               facet.by = "dose", short.panel.labs = FALSE)+
  stat_compare_means(aes(group = color),method = "wilcox.test",paired = F, label = "p.signif")+#, label = "p.signif"
  labs(x=NULL,y="proportion")
  # scale_y_continuous(limits = c(0,1.001),expand = c(0,0))+
ggsave("boxplot.png", p, width = 12, height = 7,units = "in",dpi = 600)
ggsave("boxplot.pdf", p, width = 12, height = 7,dpi = 600)

gene <- read.csv('./03.模型构建/risk_GSE98793.csv',row.names = 1)
data2 <- gene[,1:8]
data1 <- raw
