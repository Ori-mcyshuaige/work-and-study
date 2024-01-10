library(openxlsx)
library(ggplot2)
library(ggrepel)###解决标签重叠问题
library(pals)###jet颜色配色
library(stringr)
setwd("E:/DSCDY220930/LY23020703——DSCDY220930测序/3.miRNA测序分析报告/01.表格/DSCDY220930_miRNA_FGG/FGG差异miRNA")
# data<-read.xlsx("dif_group_CRC_vs_Control.xlsx",sheet = 1,colNames = T,rowNames = F,check.names = F,sep.names = " ")
data<-read.csv("./01.差异分析/genestat.csv",header = TRUE,as.is = TRUE,encoding = "UTF-8")

data$type<-"not significant"
data$type[data$pvalue<0.05&data$logfc<0]<-"down"
data$type[data$pvalue<0.05&data$logfc>0]<-"up"
data$type[data$vip<1]<-"not significant"
# data$logp<-(0-log10(data$`SLE-Control_P-VALUE`))
VolcanoPlot<-data[,c("X1","pvalue","fc","logfc","type","vip")]
colnames(VolcanoPlot)<-c("name","p","fc","logfc","type","vip")
VolcanoPlot$type<-factor(VolcanoPlot$type,levels = c('down', 'not significant', 'up'))
p<-ggplot()+
  geom_blank(data = VolcanoPlot,aes(-logfc, -log10(p)))+
  geom_hline(
    yintercept = c(-log10(0.05)),
    color = "grey50",linetype = "dashed"
  )+
  geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
  geom_point(
    data = VolcanoPlot[VolcanoPlot$type=='not significant',],
    aes(logfc,-log10(p), color = type,size = vip)
  )+
  geom_point(
    data = VolcanoPlot[VolcanoPlot$type!='not significant',],
    aes(logfc,-log10(p), color = type,size = vip)
  )+
  scale_color_manual(
    name = "Status", values = c(
      'down' = '#619cffa0', 'not significant' = '#b3b3b350', 
      'up' = '#f8766da0'
    ), guide = guide_legend(order = 1, override.aes = list(size = 3))
  )+
  scale_size_continuous(
    name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6),
    breaks = c(min(VolcanoPlot$vip), max(VolcanoPlot$vip)),
    labels = c('0.0', round(max(VolcanoPlot$vip), 1))
  )+
  theme_bw()+
  theme(
    # legend.key = element_blank()
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",size = 1),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.key.size = unit(0.8,"cm"),
    axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
  )+
  labs(
    x = expression(paste(plain("log"["2"]), ' Fold Change')), 
    y = expression(paste(plain("-log"["10"]), ' ', italic('P'), '-value'))
  )
ggsave(
  'volcano plot.jpg', p,
  width = 10, height = 7.5, units = 'in', dpi = 600
)
ggsave(
  'volcano plot.pdf', p,
  width = 10, height = 7.5, units = 'in', dpi = 600
)
# jpeg('volcano plot.jpg',width = 6000,height = 3750,res = 600)
# p
# dev.off()
# pdf('volcano plot.pdf',width = 12,height = 7.5)
# p
# dev.off()

p<-ggplot()+
  geom_blank(data = VolcanoPlot,aes(-logfc, -log10(p)))+
  geom_hline(
    yintercept = c(-log10(0.05)),
    color = "grey50",linetype = "dashed"
  )+
  geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
  geom_point(
    data = VolcanoPlot,
    aes(logfc,-log10(p), color = type),size = 3.5,alpha=0.3
  )+
  scale_color_manual(
    name = "Status", values = color, guide = guide_legend(order = 1, override.aes = list(size = 3))
  )+
  theme_bw()+
  theme(
    # legend.key = element_blank()
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",size = 1),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.key.size = unit(0.8,"cm"),
    axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
  )+
  labs(
    x = expression(paste(plain("log"["2"]), ' Fold Change')), 
    y = expression(paste(plain("-log"["10"]), ' ', italic('P'), '-value'))
  )



####无vip值
data <- read.csv('./01.差异分析/genestat.csv',header = T)
colnames(data) <- c('X1','logfc','ave','t','pvalue','adjp','b')
data$type<-"not significant(686)"
data$type[data$padj<0.05&data$log2FoldChange<(-1)]<-"down(71)"
data$type[data$padj<0.05&data$log2FoldChange>1]<-"up(48)"
# data$type[data$vip<1]<-"not significant"
# data$logp<-(0-log10(data$`SLE-Control_P-VALUE`))
VolcanoPlot<-data[,c("X","padj","log2FoldChange","type",'label')]
colnames(VolcanoPlot)<-c("name","p","logfc","type",'label')
VolcanoPlot$type<-factor(VolcanoPlot$type,levels = c('down(71)', 'not significant(686)', 'up(48)'))
VolcanoPlot$p[is.na(VolcanoPlot$p)] <- 1
p<-ggplot()+
  geom_blank(data = VolcanoPlot,aes(-logfc, -log10(p)))+
  geom_hline(
    yintercept = c(-log10(0.05)),
    color = "grey50",linetype = "dashed"
  )+
  geom_vline(xintercept = 1, color = 'grey50', linetype = 'dashed') +
  geom_vline(xintercept = -1, color = 'grey50', linetype = 'dashed') +
  geom_point(
    data = VolcanoPlot[VolcanoPlot$type=='not significant(686)',],
    aes(logfc,-log10(p), color = type),size=5
  )+
  geom_point(
    data = VolcanoPlot[VolcanoPlot$type!='not significant(686)',],
    aes(logfc,-log10(p), color = type),size=5
  )+
  scale_color_manual(
    name = 'up-down', values = c(
      'down(71)' = '#619cffa0', 
      'not significant(686)' = 'grey50',
      'up(48)' = '#f8766da0'
    ), guide = guide_legend(order = 1, override.aes = list(size = 5))
  )+
  # geom_text(data=VolcanoPlot,mapping = aes(-logfc, -log10(p),label=label))+
  geom_text_repel(data=VolcanoPlot,
                  mapping = aes(logfc, -log10(p),label=label),
                  size=5,color="black",direction="both",
                  min.segment.length = 0.05,
                  segment.alpha=0.6,label.padding = 0.4,
                  max.overlaps =100,nudge_x = 0.2,nudge_y=0.2)+
  theme_bw()+
  theme(
    # legend.key = element_blank()
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",linewidth = 1),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.key.size = unit(0.8,"cm"),
    axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
  )+
  labs(
    x = expression(paste(plain("log"["2"]), ' Fold Change')), 
    y = expression(paste(plain("-log"["10"]), ' ', italic('P'), '-value'))
  )
ggsave(
  'volcano plot.jpg', p,
  width = 10, height = 7.5, units = 'in', dpi = 600
)
ggsave(
  'volcano plot.pdf', p,
  width = 10, height = 7.5, units = 'in', dpi = 600
)
# jpeg('volcano plot.jpg',width = 6000,height = 3750,res = 600)
# p
# dev.off()
# pdf('volcano plot.pdf',width = 12,height = 7.5)
# p
# dev.off()

data <- read.csv('genestat_GSE23878.csv',header = T)
colnames(data) <- c('X1','logfc','ave','t','pvalue','adjp','b')
data$type<-"not significant"
data$type[data$adjp<0.05&data$logfc<(-log2(1.5))]<-"down"
data$type[data$adjp<0.05&data$logfc>log2(1.5)]<-"up"

# data$type[data$vip<1]<-"not significant"
# data$logp<-(0-log10(data$`SLE-Control_P-VALUE`))
VolcanoPlot<-data[,c("X1","adjp","logfc","type")]
colnames(VolcanoPlot)<-c("name","p","logfc","type")
VolcanoPlot$type<-factor(VolcanoPlot$type,levels = c('down', 'not significant', 'up'))
VolcanoPlot$p[is.na(VolcanoPlot$p)] <- 1
write.csv(VolcanoPlot,'VolcanoPlot_GSE23878.csv')
p<-ggplot()+
  geom_blank(data = VolcanoPlot,aes(-logfc, -log10(p)))+
  geom_hline(
    yintercept = c(-log10(0.05)),
    color = "grey50",linetype = "dashed"
  )+
  geom_vline(xintercept = log2(1.5), color = 'grey50', linetype = 'dashed') +
  geom_vline(xintercept = -log2(1.5), color = 'grey50', linetype = 'dashed') +
  geom_point(
    data = VolcanoPlot[VolcanoPlot$type=='not significant',],
    aes(logfc,-log10(p), color = type),size=5
  )+
  geom_point(
    data = VolcanoPlot[VolcanoPlot$type!='not significant',],
    aes(logfc,-log10(p), color = type),size=5
  )+
  scale_color_manual(
    name = 'up-down', values = c(
      'down' = '#619cffa0', 
      'not significant' = 'grey50',
      'up' = '#f8766da0'
    ), guide = guide_legend(order = 1, override.aes = list(size = 5))
  )+
  # geom_text(data=VolcanoPlot,mapping = aes(-logfc, -log10(p),label=label))+
  theme_bw()+
  theme(
    # legend.key = element_blank()
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",linewidth = 1),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.key.size = unit(0.8,"cm"),
    axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
  )+
  labs(
    x = expression(paste(plain("log"["2"]), ' Fold Change')), 
    y = expression(paste(plain("-log"["10"]), ' ', italic('adjP'), '-value'))
  )
ggsave(
  'volcano plot.jpg', p,
  width = 10, height = 7.5, units = 'in', dpi = 600
)
ggsave(
  'volcano plot.pdf', p,
  width = 10, height = 7.5, units = 'in', dpi = 600
)

  