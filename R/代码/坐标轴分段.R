library(openxlsx)
library(ggplot2)
library(ggrepel)###解决标签重叠问题
library(pals)###jet颜色配色
library(stringr)
library(ggpubr)
setwd("D:/马驰宇/181/戴主任/工作/生信流程/项目/分析/已交付AQ-LF20200701/AQ-LF20200701-QTOF-NJJM个性化散点图数据和需求(1)/AQ-LF20200701-QTOF-NJJM个性化散点图数据和需求")
data<-read.xlsx("负-lipid-筛选.xlsx",sheet = 1,colNames = T,rowNames = F,check.names = F,sep.names = " ")

data$type<-"Others"
data$type[data$`important`=="PE"]<-"PE"
data$type[data$`important`=="PC"]<-"PC"
data$type[data$`important`=="Cer-AP"]<-"Cer-AP"
data$type[data$`important`=="SQDG"]<-"SQDG"
data$type[data$`important`=="TAG"]<-"TAG"
data$type[data$`important`=="DGTS"]<-"DGTS"
data$type[data$`important`=="DAG"]<-"DAG"
data$type[data$`important`=="FA"]<-"FA"

# data$type[data$`P-VALUE`>=0.05|data$VIP<1]<-"not significant"
# data$logp<-(0-log10(data$`SLE-Control_P-VALUE`))
VolcanoPlot<-data[,c("id","MS2 name","important","P-VALUE","LOG_FOLDCHANGE","type")]
colnames(VolcanoPlot)<-c("id","MS2 name","species","p","logfc","type")
VolcanoPlot$type<-factor(VolcanoPlot$type,levels = c("Others","PE","PC","Cer-AP","SQDG","TAG","DGTS","DAG","FA"))
# VolcanoPlot1<-VolcanoPlot[VolcanoPlot$Rd<3000,]
# VolcanoPlot2<-VolcanoPlot[VolcanoPlot$Rd>3000,]
p1<-ggplot()+
  geom_blank(data = VolcanoPlot,aes(-log2(p), Rd))+
  geom_hline(
    yintercept = 0,
    color = "black"##,linetype = "dashed"
  )+
  # geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
  geom_point(
    data = VolcanoPlot[VolcanoPlot$`type`=='Others',],
    aes(-log2(p),Rd,color = type, shape = type), size = 2
  )+
  geom_point(
    data = VolcanoPlot[VolcanoPlot$`type`!='Others',],
    aes(-log2(p),Rd,color = type, shape = type), size = 2
  )+
  scale_color_manual(
    name = "Lipid species", values = c(
      "Others" = "grey", "PE" = "#90EE90", "PC" = "#8470FF", "Cer-AP" = "#00CED1", "SQDG" = "#8B658B"
    ), guide = guide_legend(order = 1, override.aes = list(size = 3))
  )+
  scale_shape_manual(
    name = "Lipid species", values = c(
      "Others" = 16, "PE" = 17, "PC" = 17, "Cer-AP" = 17, "SQDG" = 17
    ), guide = guide_legend(order = 1, override.aes = list(size = 3))
  )+
  # scale_size_continuous(
  #   name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6),
  #   breaks = c(min(VolcanoPlot$vip), max(VolcanoPlot$vip)),
  #   labels = c('0.0', round(max(VolcanoPlot$vip), 1))
  # )+
  theme_bw()+
  theme(
    # legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # panel.border = element_rect(colour = "black",size = 1),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.7, hjust = 0),
    legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.7, hjust = 0),
    legend.key.size = unit(0.8,"cm"),
    axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.line = element_line(colour='black')
  )+
  coord_cartesian(ylim = c(0,2500),xlim = c(0,41))+
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x = expression(paste(plain("-log"["2"]), ' Pvalue')),
    y = expression(paste(plain("log"["2"]), ' Fold Change'))
  )
p2<-ggplot()+
  geom_blank(data = VolcanoPlot,aes(-log2(p), Rd))+
  # geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
  geom_point(
    data = VolcanoPlot[VolcanoPlot$`type`=='Others',],
    aes(-log2(p),Rd,color = type, shape = type), size = 2
  )+
  geom_point(
    data = VolcanoPlot[VolcanoPlot$`type`!='Others',],
    aes(-log2(p),Rd,color = type, shape = type), size = 2
  )+
  scale_color_manual(
    name = "Lipid species", values = c(
      "Others" = "grey", "PE" = "#90EE90", "PC" = "#8470FF", "Cer-AP" = "#00CED1", "SQDG" = "#8B658B"
    ), guide = guide_legend(order = 1, override.aes = list(size = 3))
  )+
  scale_shape_manual(
    name = "Lipid species", values = c(
      "Others" = 16, "PE" = 17, "PC" = 17, "Cer-AP" = 17, "SQDG" = 17
    ), guide = guide_legend(order = 1, override.aes = list(size = 3))
  )+
  # scale_size_continuous(
  #   name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6),
  #   breaks = c(min(VolcanoPlot$vip), max(VolcanoPlot$vip)),
  #   labels = c('0.0', round(max(VolcanoPlot$vip), 1))
  # )+
  theme_bw()+
  theme(
    # legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    # panel.border = element_rect(colour = "black",size = 1),
    axis.line.y = element_line(colour='black'),
    axis.text.y = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.7, hjust = 0),
    legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.7, hjust = 0),
    legend.key.size = unit(0.8,"cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  coord_cartesian(ylim = c(5000,33621.9),xlim = c(0,41))+
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x = NULL,
    y = NULL
  )
p<-ggarrange(p2,p1,heights=c(1/5, 4/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")
# jpeg('- M1 Vs Wild.jpg',width = 6000,height = 3750,res = 600)
# p
# dev.off()
# pdf('- M1 Vs Wild.pdf',width = 12,height = 7.5)
# p
# dev.off()
# ggsave(
#   'volcano plot.jpg', p, 
#   width = 12, height = 7.5, units = 'in', dpi = 600
# )
# ggsave(
#   'volcano plot.pdf', p, 
#   width = 12, height = 7.5, units = 'in', dpi = 600
# )



# p<-ggplot()+
#   geom_blank(data = VolcanoPlot,aes(-log2(p), logfc))+
#   geom_hline(
#     yintercept = 0,
#     color = "black",##,linetype = "dashed"
#     size = 1
#   )+
#   # geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
#   geom_point(
#     data = VolcanoPlot[VolcanoPlot$`type`=='Others',],
#     aes(-log2(p),logfc,color = type, shape = type), size = 3.5
#   )+
#   geom_point(
#     data = VolcanoPlot[VolcanoPlot$`type`!='Others',],
#     aes(-log2(p),logfc,color = type, shape = type), size = 3.5
#   )+
#   scale_color_manual(
#     name = "Lipid species", values = c(
#       "Others" = "grey", "PE" = "#90EE90", "PC" = "#8470FF", "Cer-AP" = "#00CED1", "SQDG" = "#8B658B", "TAG" = "#FFD700", "DGTS" = "#BC8F8F", "DAG" = "#D2691E", "FA" = "#FF4040"
#     ), guide = guide_legend(order = 1, override.aes = list(size = 3))
#   )+
#   scale_shape_manual(
#     name = "Lipid species", values = c(
#       "Others" = 16, "PE" = 17, "PC" = 17, "Cer-AP" = 17, "SQDG" = 17, "TAG" = 17, "DGTS" = 17, "DAG" = 17, "FA" = 17
#     ), guide = guide_legend(order = 1, override.aes = list(size = 3))
#   )+
#   # scale_size_continuous(
#   #   name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6),
#   #   breaks = c(min(VolcanoPlot$vip), max(VolcanoPlot$vip)),
#   #   labels = c('0.0', round(max(VolcanoPlot$vip), 1))
#   # )+
#   theme_bw()+
#   theme(
#     # legend.key = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     # panel.border = element_rect(colour = "black",size = 1),
#     legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.7, hjust = 0),
#     legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.7, hjust = 0),
#     legend.key.size = unit(0.8,"cm"),
#     axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
#     axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
#     axis.line = element_line(colour='black',size = 1),
#     axis.ticks = element_line(size = 1),
#     axis.ticks.length = unit(0.2,"cm")
#   )+
#   coord_cartesian(xlim = c(0,30.5))+
#   scale_x_continuous(expand = c(0, 0))+
#   labs(
#     x = expression(paste(plain("-log"["2"]), ' Pvalue')),
#     y = expression(paste(plain("log"["2"]), ' Fold Change'))
#   )
# ggsave('MTPf0-4 Vs Wild.jpg',p,width = 12,height = 7.5,units = 'in',dpi = 600)
# ggsave('MTPf0-4 Vs Wild.pdf',p,width = 12,height = 7.5,units = 'in',dpi = 600)
# jpeg('M1 Vs Wild.jpg',width = 6000,height = 3750,res = 600)
# p
# dev.off()
# pdf('M1 Vs Wild.pdf',width = 12,height = 7.5)
# p
# dev.off()
