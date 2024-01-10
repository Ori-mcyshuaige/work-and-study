library(dior)
scRNA = dior::read_h5(file='tcelldata.h5', target.object = 'seurat')
#####设置对比关系
group=levels(factor(scRNA$tissue_type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

library(ggplot2)
library(ggpubr)  # Assuming you are using ggpubr for theme_pubr function

desired_Celltype <- scRNA$Cell.types
filtered_scRNA <- scRNA@meta.data %>%
  filter(Cell.types %in% desired_Celltype)

p <- ggplot(filtered_scRNA, aes(x = tissue_type, y = FXN, fill = tissue_type)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, hjust = 1),
        axis.line.x = element_line(size = 0.8, linetype = "solid"),
        axis.line.y = element_line(size = 0.8, linetype = "solid"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~Cell.types) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1.5) +
  # geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     size = 4,
                     label.y = NULL,
                     hide.ns = T) +
  scale_fill_manual(values = c("#1F78B4", '#33CC00')) +##, "#FF7F00", '#FF0033'
  rotate_x_text(angle = 45) +
  ylab("FXN") +
  guides(fill = "none") +
  expand_limits(y = 0.3) 
p
ggsave('./figure7/scRNA_tissue_FXN.pdf',p,width = 9,height = 9,dpi = 600)


#####设置对比关系
group=levels(factor(scRNA$copykat.pred))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

library(ggplot2)
library(ggpubr)  # Assuming you are using ggpubr for theme_pubr function

desired_Celltype <- scRNA$Cell.types
filtered_scRNA <- scRNA@meta.data %>%
  filter(Cell.types %in% desired_Celltype)

p <- ggplot(filtered_scRNA, aes(x = copykat.pred, y = disulfidptosis_score, fill = copykat.pred)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, hjust = 1),
        axis.line.x = element_line(size = 0.8, linetype = "solid"),
        axis.line.y = element_line(size = 0.8, linetype = "solid"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~Cell.types) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1.5) +
  # geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     size = 4,
                     label.y = 0.6,
                     hide.ns = T) +
  scale_fill_manual(values = c("#1F78B4", '#33CC00')) +##, "#FF7F00", '#FF0033'
  rotate_x_text(angle = 45) +
  ylab("Disulfidptosis Score") +
  guides(fill = "none") +
  expand_limits(y = 0.3) +
  ylim(0, 0.8)
p
ggsave('./08.celltype_diff/scRNA_copykat.pdf',p,width = 9,height = 9,dpi = 600)


#####设置对比关系
group=levels(factor(scRNA$group))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

library(ggplot2)
library(ggpubr)  # Assuming you are using ggpubr for theme_pubr function

desired_Celltype <- scRNA$Celltype
Celltype <- scRNA$Celltype
filtered_scRNA <- scRNA@meta.data %>%
  filter(Celltype %in% desired_Celltype)
filtered_scRNA$Nr3c1 <- scRNA[['RNA']]@data['Nr3c1',][rownames(filtered_scRNA)]
filtered_scRNA$Rack1 <- scRNA[['RNA']]@data['Rack1',][rownames(filtered_scRNA)]
p <- ggplot(filtered_scRNA, aes(x = group, y = Nr3c1, fill = group)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, hjust = 1),
        axis.line.x = element_line(size = 0.8, linetype = "solid"),
        axis.line.y = element_line(size = 0.8, linetype = "solid"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~Celltype) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1) +
  # geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  stat_compare_means(comparisons = my_comparisons,method = 't.test',
                     
                     size = 4,
                     label.y = 2.5,
                     hide.ns = T) +
  scale_fill_manual(values = c("#1F78B4", '#33CC00')) +##, "#FF7F00", '#FF0033'
  rotate_x_text(angle = 45) +
  ylab("Nr3c1") +
  guides(fill = "none") +
  expand_limits(y = 0.3) +
  ylim(0, 3)
p
ggsave('./Nr3c1_group.pdf',p,width = 9,height = 9,dpi = 600)

# method默认为“wilcox.test”（非参数检验），可指定method = “t.test”，表示T检验（参数检验）
# 
# 返回值为具有以下列的数据框：
# 
# .y.：用于统计检验的数值变量
# 
# p：P值
# 
# p.adj：调整后的P值，调整P值的默认方法为p.adjust.method = “holm”
# 
# p.format: 格式化的P值
# 
# p.signif：显著性水平，即用不同数量的 * 表示显著性水平
# 
# method：用于组间比较的统计方法