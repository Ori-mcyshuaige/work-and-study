library(magrittr)
library(plyr)
library(openxlsx)
library(GSVA)
library(epitools) # 用于生成期望频数表
library(reshape2)
library(ggplot2)

meta <- scRNA@meta.data

observe.data <- as.matrix(as.data.frame.matrix(table(meta$group, meta$Cell.types)))
expected.data <- expected(observe.data)
Ratio <- observe.data/expected.data
plot.data <- melt(Ratio)
colnames(plot.data) <- c("Group", "Celltype", "Ratio")
plot.data$Group <- factor(plot.data$Group, levels = c("cirrhotic", "HCC"))
color <- setNames(object = c("#1F78B4","#E31A1C"), 
                  nm = c("cirrhotic", "HCC"))
p2 <-  ggplot(plot.data, aes(x = Celltype, y = Ratio, group = Group, color = Group)) +
  geom_point(size = 3) +
  geom_line(linetype = "dashed") +
  scale_color_manual(values = color) +
  theme_classic() +
  labs(x = "", y = "Ro/e", fill = "position") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, color = "black"),
    axis.title.y = element_text(hjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(hjust = 0.5, size = 16, color = "black"),
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(.9, 0.9),  # 调整图例位置到右侧中间
    legend.text = element_text(size = 14),  # 调整图例字体大小为16
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.margin = margin(t = 0, r = 4, b = 0, l = 0, unit = "lines")  # 设置右边边距为4个单位
  )
ggsave('./07.RoeDot/RoeDot.pdf',width = 9,height = 9,p2,dpi = 600)
p2
