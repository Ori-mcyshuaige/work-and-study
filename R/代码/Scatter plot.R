library(ggpubr)# 载入ggpubr这个包
library(openxlsx)
library(pals)

# setwd("C:\\Users\\dujie\\Desktop") #别忘记进入工作目录，就是你放文件的地方
data <-read.xlsx("spearman data.xlsx",sheet = 2,rowNames = T)
a <- names(data)
gene <- a[1]
meta <- a[2]
for(i in gene){
  for (j in meta){
    plotname <- paste(i, "-", j, sep = "")
    print(plotname)
    p <- ggscatter(data, x = i,y = j,shape = 20,size = 5,
                                 color = "#FF0000",cor.method = "spearman",
                                 add.params = list(color = "black"),
                                 add = "reg.line",#label = "symbol",label.select = data$symbol,
                                 conf.int = T, cor.coef = T,show.legend.text = TRUE)+
      theme_bw() +
      theme(
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
    # corr <- cor.test(data[,i],data[,j],method ="spearman")$estimate
    # Pvalue <- cor.test(data[,i],data[,j],method ="spearman")$p.value
    #  p <- p + annotate("text",x = max(data[,i])*0.8, y = max(data[,j])*1.05,label=paste0('rho = ', round(corr,3))) +
    #   annotate("text",x = max(data[,i])*0.8, y = max(data[,j])*1.02,label=paste0('p = ', round(Pvalue,3)))   
	  ggsave(paste(plotname,".jpg", sep = ""), p, width = 6, height = 7.5, units = 'in', dpi = 600)
	  ggsave(paste(plotname,".pdf", sep = ""),p,width = 6,height = 7.5,units = 'in',dpi = 600)
  }
 }
 
raw<-read.xlsx("SLEDAI-挑选蛋白-Corrected数据-散点图.xlsx",sheet=2)
data<-data.frame()
for(i in colnames(raw)[2:9]){
  line<-data.frame(raw[,1],raw[,i])
  data<-rbind(data,line)
}
name<-sapply(colnames(raw)[2:9],function(i){
  strsplit(i,"-")[[1]][1]
}
)
Group<-rep(name,each=121)
data$Protein<-Group
colnames(data)[1:2]<-c("x","y")
color<-jet(8)
names(color)<-unique(data$Protein)
p<-ggplot(data=data,mapping=aes(x=x,y=y,colour=Protein))+
  geom_point()+
  scale_color_manual(values = color)+
  stat_smooth(method = "lm",alpha=0.2)+
  theme_bw() +
  theme(
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
  )+
  labs(x="SLEDAI-adjusted",y="Protein-adjusted")
ggsave("偏相关散点.png",p, width = 4.8, height = 4.8,units = "in",dpi = 600)
ggsave("偏相关散点.pdf",p, width = 4.8, height = 4.8)
