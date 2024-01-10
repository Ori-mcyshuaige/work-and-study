library(ggplot2)
library(ggrepel)###解决标签重叠问题
library(pals)###jet颜色配色
library(stringr)
library(openxlsx)
 


#####自定义气泡图
ggplot() +
  geom_blank(data = data, aes(Impact, `-ln(p)`, size = Impact)) +
  geom_point(
    data = data[data$`Raw p` < 0.1,],
    aes(Impact, `-ln(p)`, size = Impact),col = "blue"
  )+
  geom_point(
    data = data[data$`Raw p` >= 0.1,],
    aes(Impact, `-ln(p)`, size = Impact),col = "blue",alpha = 0.4
  )+
  geom_text_repel(
    data = data,
    mapping = aes(Impact, `-ln(p)`,  label = Pathway,color=p),
    size = 5,family = "serif",color = "black",
    point.padding = unit(1, 'lines'),
    min.segment.length = unit(0.1, "lines"),
    segment.size = 0.6,
    # segment.color = 'green',
    force=7,box.padding=0.2,
    show.legend = F
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    legend.text = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(color = 'black', size = 1.5),
    panel.grid.minor = element_blank(),
    # legend.text = element_text(size = 16, color = "black",family = "serif",
    #                            face = "bold", vjust = 0.5, hjust = 0),
    # legend.title = element_text(size = 16, color = "black",family = "serif",
    #                             face = "bold", vjust = 0.5, hjust = 0),
    # legend.key.size = unit(0.8,"cm"),
    axis.text = element_text(size = 16, color = "black",family = "serif",
                             face = "bold", vjust = 0.5, hjust = 0.5),  ###  控制坐标轴刻度的大小颜色
    axis.title = element_text(size = 16, color = "black", family = "serif",
                              face = "bold", vjust = 0.5, hjust = 0.5) ###  控制坐标轴标题的大小颜色
  ) +
  scale_size_continuous(
    name = 'Impact', range = c(8, 16),
    breaks = c(min(data$Impact), max(data$Impact)),
    labels = c(min(data$Impact), round(max(data$Impact), 3)),
    guide = guide_legend(title.hjust = 0.5)
  )+
  labs(
    x = 'Impact',
    y = expression(paste(-ln, ' ', italic('P'), '-value'))
  )

####百趣气泡图
setwd("D:/马驰宇/181/戴主任/工作/生信流程/项目/分析/BQ-YP20191115-FC4-2")
df.ma<-read.csv("pathway_results2.csv",row.names = NULL)
df.ma<-df.ma[,c(1,4,5,8,9)]
colnames(df.ma)<- c('pathway','hits','p','fdr','impact')
df.ma$pathway<-factor(df.ma$pathway,levels = df.ma$pathway)
# if (length(unique(df.ma$impact))<=1) {
#   next()
# }
if(length(which(df.ma$impact > 0)) >= 5) {
  filter <- expression(
    df.ma$impact >= sort(df.ma$impact, decreasing = T)[5]
  )
} else {
  filter <- expression(df.ma$impact >= 0)
}
p <- ggplot() +
  geom_point(
    data = df.ma, shape = 21, color = 'black', stroke = 1, 
    mapping = aes(x = impact, y = -log(p,10), size = impact, fill = p)
  ) +
  geom_label_repel(
    data = df.ma[eval(filter) | df.ma$`p` < 0.05, ], 
    mapping = aes(impact, -log(p,10), fill = p, label = pathway), 
    color = 'black', size = 3,
    label.padding = unit(0.2, 'lines'), 
    point.padding = unit(0.5, 'lines'), 
    min.segment.length = unit(0.1, "lines"), 
    segment.color = 'grey50', segment.size = 1, 
    show.legend = F
  ) +
  geom_hline(
    yintercept = c(-log10(0.05)),
    color = "grey50",linetype = "dashed"
  )+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",size = 1),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.key.size = unit(0.8,"cm"),
    axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
  ) +
  scale_size_continuous(
    name = 'Impact', range = c(8, 16), 
    breaks = c(min(df.ma$impact), max(df.ma$impact)), 
    labels = c(min(df.ma$impact), round(max(df.ma$impact), 3)), 
    guide = guide_legend(title.hjust = 0.5)
  ) + 
  scale_fill_gradient(
    name = "P-value", 
    low = '#FF4040', high = '#86ABED', 
    guide = guide_colorbar(
      title.position = 'top', title.hjust = 0.5, 
      direction = 'horizontal'
    )
  ) +
  labs(
    x = 'Impact', 
    y = expression(paste(-log10, ' ', italic('P'), '-value'))
  )
ggsave(
  'Bubble Plot2.jpg', p,
  width = 7.2, height = 5.4, units = 'in', dpi = 600
)
ggsave(
  'Bubble Plot2.pdf', p,
  width = 7.2, height = 5.4
)





###常规富集分析气泡图
setwd("C:/Users/Mcy-aiyu/Desktop/何院专辑/研究院项目/戴主任项目重分析/代谢分析/Go")
data<-read.xlsx("analysis-GO-MF(1).xlsx",sheet = 3,colNames = T,rowNames = F,check.names = F,sep.names = " ")

x=data$`upload_1 (fold Enrichment)`
y=factor(data$`GO molecular function complete`,levels = data$`GO molecular function complete`)

p<-ggplot(data,aes(x,y))+
  geom_point(aes(size = `upload_1 (88)`,color = -log10(`upload_1 (FDR)`))) +
  scale_color_gradientn(values = c(0,1),colors = c("#4876FF","red1")) +
  scale_size_continuous(
    name = 'Hits Number', range = c(4, 12),
    breaks = c(min(data$`upload_1 (88)`), max(data$`upload_1 (88)`)),
    labels = c(min(data$`upload_1 (88)`), round(max(data$`upload_1 (88)`), 4)),
    guide = guide_legend(title.hjust = 0.5)
  )+
  theme_bw()+
  theme(
    panel.grid.major = element_line(
      color = 'grey', linetype = 'dashed'
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )+
  labs(
    x = "Fold Enrichment",
    y = "GO molecular function complete",
    size = "Hits Number",
    color = "-log10(FDR)"
  )

p<-ggplot(data,aes(x,y))+
  geom_point(aes(size = `upload_1 (88)`,color = color)) +
  scale_color_manual(values = c("<0.05"="#FFB55F","<0.01"="#FF7735","<0.001"="#FF0005")) +

  scale_size_continuous(
    name = 'Hits Number', range = c(4,12),
    breaks = c(min(data$`upload_1 (88)`), max(data$`upload_1 (88)`)),
    labels = c(min(data$`upload_1 (88)`), round(max(data$`upload_1 (88)`), 4)),
    guide = guide_legend(title.hjust = 0.5)
  )+
  theme_bw()+
  theme(
    panel.grid.major = element_line(
      color = 'grey', linetype = 'dashed'
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    # legend.key.size = unit(2,"cm")
  )+
  labs(
    x = "Fold Enrichment",
    y = "GO molecular function complete",
    size = "Hits Number",
    color = "FDR"
  )


ggsave(
  'Bubble Plot.jpg', p,
  width = 11.5, height = 10, units = 'in', dpi = 600
)
ggsave(
  'Bubble Plot.pdf', p,
  width = 11.5, height = 10
)

  
##又自定义了
p<-ggplot() +
  geom_point(
    data = data,stroke=1,shape=21,
    aes(Impact, `-ln(p)`, size = Impact,fill=`-ln(p)`)
  )+
  geom_label_repel(
    data = data,
    mapping = aes(Impact, `-ln(p)`, label=Pathway,fill=`-ln(p)`),
    size = 5,family = "sans",color = "black",
    point.padding = unit(1, 'lines'),
    min.segment.length = unit(0.1, "lines"),
    segment.size = 0.6,
    # segment.color = 'green',
    force=7,box.padding=0.2,
    show.legend = F
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(
      color = '#3f00f6', linetype = 'dashed'
    ), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'black', size = 1.5),
    legend.text = element_text(size = 9, color = "black",family = "sans",
                               vjust = 0.5, hjust = 0.8),
    legend.title = element_text(size = 14, color = "black",family = "sans",
                                vjust = 0.5, hjust = 0),
    # legend.key.size = unit(0.8,"cm"),
    axis.text = element_text(size = 14, color = "black",family = "sans",
                             vjust = 0.5, hjust = 0.5),  ###  控制坐标轴刻度的大小颜色
    axis.title = element_text(size = 14, color = "black", family = "sans",
                              vjust = 0.5, hjust = 0.5) ###  控制坐标轴标题的大小颜色
  ) +
  scale_size_continuous(
    name = 'Impact', range = c(8, 16),
    breaks = c(min(data$logFC), max(data$logFC)),
    labels = c(min(data$logFC), round(max(data$logFC), 3)),
    guide = guide_legend(title.hjust = 0.5)
  )+
  scale_fill_gradient(
    name = "-ln P-value", 
    low = '#86ABED', high = '#FF4040', 
    guide = guide_colorbar(
      title.position = 'top', title.hjust = 0.5, 
      direction = 'horizontal'
    )
  ) +labs(
    x = 'Impact',
    y = "-ln P-value"
  )

ggsave(
  'Bubble Plot.jpg', p,
  width = 9.6, height = 7.2, units = 'in', dpi = 600
)
ggsave(
  'Bubble Plot.pdf', p,
  width = 9.6, height = 7.2
)


###第N个气泡图
pathwayplot<-gsearesult[pnumber,]
pathwayplot<-pathwayplot[order(pathwayplot$NES),]
pathwayplot$ID<-factor(pathwayplot$ID,levels = pathwayplot$ID)
color<-c("#FF3030","#FF303050","#4876FF","#4876FF50","#696969")
names(color)<-c("activated(p.adjust<0.05)","activated(p<0.05)","suppressed(p.adjust<0.05)","suppressed(p<0.05)","not significant")

p<-ggplot(pathwayplot,aes(NES,ID))+
  geom_point(
    aes(NES, ID, size = count,col = color)
  )+
  scale_color_manual(values = color)+
  scale_size_continuous(
    name = 'Hits Number', range = c(4, 12),
    breaks = c(min(pathwayplot$count), max(pathwayplot$count)),
    labels = c(min(pathwayplot$count), round(max(pathwayplot$count), 4)),
    guide = guide_legend(title.hjust = 0.5)
  )+
  theme_bw()+
  theme(
    panel.grid.major = element_line(
      color = 'grey', linetype = 'dashed'
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.position = c(.95, .05),
    legend.justification = c(.95, .05)
  )+
  labs(
    x = "NES",
    y = "Pathway",
    size = "Hits Number"
  )
ggsave(paste0(temp1,"/",temp2,"/","pathwayplot.png"),p,dpi = 600,width = 0.75*nrow(pathwayplot),height = 9,limitsize = F,units = "in")
ggsave(paste0(temp1,"/",temp2,"/","pathwayplot.pdf"),p,width = 0.75*nrow(pathwayplot),height = 9,limitsize = F,units = "in")
