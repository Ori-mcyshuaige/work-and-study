ggplot(dat, aes(x=expr, colour=Name)) + geom_density(size=1)+
theme_bw() +
  theme(
    # legend.position = "none",
    panel.grid.major = element_blank(),
    # legend.text = element_blank(),
    # legend.title = element_blank(),
    panel.border = element_rect(color = 'black', size = 1.5),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, color = "black",family = "sans",
                               face = "bold", vjust = 0.5, hjust = 0),
    legend.title = element_blank(),
    # legend.key.size = unit(0.8,"cm"),
    axis.text = element_text(size = 12, color = "black",family = "sans",
                             face = "bold", vjust = 0.5, hjust = 0.5),  ###  控制坐标轴刻度的大小颜色
    axis.title = element_text(size = 12, color = "black", family = "sans",
                              face = "bold", vjust = 0.5, hjust = 0.5) ###  控制坐标轴标题的大小颜色
  )