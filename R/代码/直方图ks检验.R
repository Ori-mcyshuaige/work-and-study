library(ggplot2)

p<-ggplot(data=allamp,aes(x=log2Arel,fill=group))+
  geom_histogram(binwidth=1,alpha=0.6,stat="count")
set.seed(135)
df <- data.frame(GFP=c(rnorm(60,5,1.5),rnorm(90,12,2)),
                 MSR=c(rnorm(110,5,1.5),rnorm(40,12,2.5)))
data <- melt(df)
p <- ggplot(allamp,aes(x=log2Arel,fill=group))+
  geom_histogram(binwidth=.15,position = 'identity',alpha=0.6
                 )+ #使用aes(y=after_stat(count/sum(count)))则可以绘制频率分布
  scale_fill_manual(values = c('red','black'))+
  scale_color_manual(values = c('red','black'))+
  geom_vline(xintercept = mean(allamp[allamp$group=='noncircadian',]$log2Arel),color = "black",linetype = "dashed")+
  geom_vline(xintercept = mean(allamp[allamp$group=='circadian',]$log2Arel),color = "red",linetype = "dashed")+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10,colour = "black"),
    axis.title = element_text(size = 10,colour = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 1,colour = "black")
  )
ggsave("histogram.png", p, width = 4.8, height = 4.8,units = "in",dpi = 600)
ggsave("histogram.pdf", p, width = 4.8, height = 4.8,dpi = 600)

mean(allamp[allamp$group=='circadian',]$log2Arel)
mean(allamp[allamp$group=='noncircadian',]$log2Arel)

ks.test(allamp[allamp$group=='circadian',]$log2Arel,allamp[allamp$group=='noncircadian',]$log2Arel)


p <- ggplot(allamp,aes(x=log2Arel,color=group))+
  geom_density()+ #使用aes(y=after_stat(count/sum(count)))则可以绘制频率分布
  # scale_fill_manual(values = c('red','black'))+
  scale_color_manual(values = c('red','black'))+
  # geom_vline(xintercept = mean(allamp[allamp$group=='noncircadian',]$log2Arel),color = "black",linetype = "dashed")+
  # geom_vline(xintercept = mean(allamp[allamp$group=='circadian',]$log2Arel),color = "red",linetype = "dashed")+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10,colour = "black"),
    axis.title = element_text(size = 10,colour = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 1,colour = "black")
  )
ggsave("density.png", p, width = 4.8, height = 4.8,units = "in",dpi = 600)
ggsave("density.pdf", p, width = 4.8, height = 4.8,dpi = 600)

p <- ggplot(allamp,aes(x=log2Arel,color=group))+
  stat_ecdf()+ #使用aes(y=after_stat(count/sum(count)))则可以绘制频率分布
  # scale_fill_manual(values = c('red','black'))+
  scale_color_manual(values = c('red','black'))+
  # geom_vline(xintercept = mean(allamp[allamp$group=='noncircadian',]$log2Arel),color = "black",linetype = "dashed")+
  # geom_vline(xintercept = mean(allamp[allamp$group=='circadian',]$log2Arel),color = "red",linetype = "dashed")+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10,colour = "black"),
    axis.title = element_text(size = 10,colour = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 1,colour = "black")
  )
ggsave("ECDF.png", p, width = 4.8, height = 4.8,units = "in",dpi = 600)
ggsave("ECDF.pdf", p, width = 4.8, height = 4.8,dpi = 600)
