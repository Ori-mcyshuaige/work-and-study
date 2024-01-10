sing_plot <- function(umap,copykat){
library(ggrepel)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
cell_type_med <- umap %>%
  group_by(copykat) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = copykat)) +
            geom_point(size = 1 , alpha =1 )  +
            scale_color_manual(values = allcolour)+
            theme(panel.grid.major = element_blank(), #主网格线
            panel.grid.minor = element_blank(), #次网格线
            panel.border = element_blank(), #边框
            axis.title = element_blank(),  #轴标题
            axis.text = element_blank(), # 文本
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = 'white'), #背景色
            plot.background=element_rect(fill="white"))+ theme(
            legend.title = element_blank(), #去掉legend.title 
            legend.key=element_rect(fill='white'), #
            legend.text = element_text(size=20), #设置legend标签的大小
            legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
            guides(color = guide_legend(override.aes = list(size=5))) + #设置legend中 点的大小  
            geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
            geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
           annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, 
                    label = "UMAP_1", color="black",size = 5, fontface="bold" ) + 
           annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, 
                    label = "UMAP_2",color="black",size = 5, fontface="bold" ,angle=90) +
           geom_label_repel(aes(label=copykat), fontface="bold",data = cell_type_med,
                           point.padding=unit(0.5, "lines")) +
           theme(legend.position = "none")
}

