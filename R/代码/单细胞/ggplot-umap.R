

library(ggplot2)
library(dplyr)
library(tidyverse)
#原文颜色
nature_col = c("#AEC7E8","#FFBB78","#9467BD","#7200DA",
               "#17BECF","#FF7F0E","#C49C94","#2CA02C","#8C564B",
               "#E377C2","#D62728","#FF9896","#98DF8A","#BCBD22","#C5B0D5")

df <- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","cluster"))

###指定显示顺序
cluster_order <- c('YSMP','ErP',"MkP","GMP",
                   "Myeloblast","Monocyte","Mac_1",
                   "Mac_2","Mac_3","Mac_4",
                   "HSPC","CD7loP","CD7hiP",'ILC',"Mast cell")
scales::show_col(nature_col)
###指定因子及其颜色、一一对应
col_cluster <- setNames(c("#AEC7E8","#FFBB78","#9467BD","#7200DA",
                          "#17BECF","#FF7F0E","#C49C94","#2CA02C","#8C564B",
                          "#E377C2","#D62728","#FF9896","#98DF8A","#BCBD22","#C5B0D5"),
                    c('YSMP','ErP',"MkP","GMP",
                      "Myeloblast","Monocyte","Mac_1",
                      "Mac_2","Mac_3","Mac_4",
                      "HSPC","CD7loP","CD7hiP",'ILC',"Mast cell"))

#作图-------UMAP设置
ggplot(df,aes(x= UMAP1 , y = UMAP2 ,col=factor(cluster, levels = cluster_order))) + 
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.8,'cm'),
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))

#作图-------UMAP坐标轴设置

ggplot(df,aes(x= UMAP1 , y = UMAP2 ,col=factor(cluster, levels = cluster_order))) + 
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP1) , y = min(UMAP2),xend = min(UMAP1)+3, yend = min(UMAP2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"), 
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP1), y = min(UMAP2),xend = min(UMAP1),yend = min(UMAP2)+1.5),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), 
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP1) +1.4, y = min(df$UMAP2) -0.3, label = "UMAP1",
           color="black",size = 5) + 
  annotate("text", x = min(df$UMAP1) -0.8, y = min(df$UMAP2) + 0.8, label = "UMAP2",
           color="black",size = 5,angle=90) 



####===========添加分组标记
#添加cluster名称
cell <- df %>%group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1),
            UMAP2 = median(UMAP2))

rownames(cell) <- cell$cluster
A <- cell[cluster_order,]
a <- c(1:15)
A$ID <- a
A$ID <- as.factor(A$ID)

p <- ggplot(df,aes(x= UMAP1 , y = UMAP2 ,col=factor(cluster, levels = cluster_order))) + 
  geom_point(size = 2, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none',
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
  



#制作legend
B <- A
B$x <- 1
B$lab <- c(15:1)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>% 
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=cluster)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 6)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)+
  annotate("text", x = 1, y = 16, label = "Cluster",
           color="black",size = 6)+
  geom_text(aes(label = cluster, x = x+0.001, 
                y = Response_round), size = 5, hjust=0)+
  scale_x_continuous(expand=c(-0.01,0.01))

#组合
library(cowplot)
plotlist <- list(p, leg)
plot_grid(plotlist = plotlist, ncol = 2, align="hv", rel_widths = c(5,2))



##########=============================================================



