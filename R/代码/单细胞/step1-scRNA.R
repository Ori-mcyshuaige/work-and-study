rm(list = ls())
options(stringsAsFactors = F)
gc()
library(tidyverse)
library(patchwork)
library(dior)
library(scater)
library(Seurat)
library(cowplot)
library(reticulate)
library(dior)
library(RColorBrewer)
library(ClusterGVis)
library(org.Hs.eg.db)
library(viridis)
library(ggpubr)
library(org.Hs.eg.db)

# col <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
#           brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
scRNA = dior::read_h5(file='./step1-scRNA/data.h5', target.object = 'seurat')
cols <- c("#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#A6CEE3", "#1F78B4", "#B2DF8A" ,"#33A02C")
mydata<- FetchData(scRNA,vars = c("UMAP_1","UMAP_2","Celltype"))
mydata$Celltype<- factor(mydata$Celltype, levels = c("CD4 T cells", "CD8 T cells","B cells",
                                                     "Macrophage","Mast cells","Epithelial cells"))
p1 <- ggplot(mydata,aes(x= UMAP_1 , y = UMAP_2 ,color = Celltype)) +geom_point(size = 0.1, alpha =0.7)+
  scale_color_manual(values = cols)+
  guides(colour=guide_legend(override.aes=list(size=5)))+ #改标签大小
  theme(axis.text.y = element_blank(),   #去掉y轴刻度注释
        axis.ticks.y = element_blank(),    #去掉y轴刻度
        axis.text.x = element_blank(),   #去掉x轴刻度注释
        axis.ticks.x = element_blank())+  #去掉x轴刻度
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+  #加边框
  theme(panel.background = element_rect(fill = "white"))+
  labs(title="Cell Type")+theme(plot.title = element_text(hjust = 0.5))
ggsave('./step1-scRNA/celltype.pdf',p1,width = 8,height = 7,dpi = 600) 

match_celltype_levels <- c("CD4 T cells", "CD8 T cells","B cells",
                           "Macrophage","Mast cells","Epithelial cells")
mydata <- scRNA@meta.data
mydata1 <- mydata[mydata$sample=='Ag',]
mydata2 <-mydata1%>%
  group_by(phenotype) %>%
  mutate(celltype = factor(Celltype,levels = match_celltype_levels)) %>%
  arrange(celltype)
mydata2$phenotype <- factor(mydata2$phenotype, levels = c("ANA", "AA"))
p2 <- ggplot() +
  geom_bar(data = mydata2, aes(x = phenotype, fill = factor(celltype)), position = position_fill(reverse = TRUE)) +scale_fill_manual(values = cols) +
  labs(fill = "Celltype",x ='Group', y = "Fraction of cells",title = 'Allergen')+
  theme(
    plot.title = element_text(hjust = 0.5), # 标题居中
    panel.grid.major = element_blank(), # 主网格线
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"), # 背景色
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")
  )

mydata3 <- mydata[mydata$sample=='Dil',]
mydata3 <-mydata3%>%
  group_by(phenotype) %>%
  mutate(celltype = factor(Celltype,levels = match_celltype_levels)) %>%
  arrange(celltype)
mydata3$phenotype <- factor(mydata3$phenotype, levels = c("ANA", "AA"))
p3 <- ggplot() +
  geom_bar(data = mydata3, aes(x = phenotype, fill = factor(celltype)), position = position_fill(reverse = TRUE)) +scale_fill_manual(values = cols) +
  labs(fill = "Celltype",x ='Group', y = "Fraction of cells",title = 'Diluent')+
  theme(
    plot.title = element_text(hjust = 0.5), # 标题居中
    panel.grid.major = element_blank(), # 主网格线
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"), # 背景色
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")
  )
mydata4 <- mydata[mydata$sample=='Pre',]
mydata4 <-mydata4%>%
  group_by(phenotype) %>%
  mutate(celltype = factor(Celltype,levels = match_celltype_levels)) %>%
  arrange(celltype)
mydata4$phenotype <- factor(mydata4$phenotype, levels = c("ANA", "AA"))
p4 <- ggplot() +
  geom_bar(data = mydata4, aes(x = phenotype, fill = factor(celltype)), position = position_fill(reverse = TRUE)) +scale_fill_manual(values = cols) +
  labs(fill = "Celltype",x ='Group', y = "Fraction of cells",title = 'Baseline')+
  theme(
    plot.title = element_text(hjust = 0.5), # 标题居中
    panel.grid.major = element_blank(), # 主网格线
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"), # 背景色
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")
  )
p5 <- p2+p3+p4+ plot_layout(guides='collect',widths = c(5, 5, 5))
ggsave('./step1-scRNA/Fraction.pdf',p5,height =8,width = 12,dpi = 500 )


#####compare
cu <- read.table('./step1-scRNA/gene.txt',sep = '\t')
gene1 <- as.list(cu)
scRNA <- AddModuleScore(scRNA,features = gene1,name = 'cuproptosis')
scRNA$phenotype <- factor(scRNA$phenotype, levels = c("ANA", "AA"))
group <- levels(factor(scRNA$phenotype))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

p6 <- scRNA@meta.data[scRNA@meta.data$sample=='Ag',] %>%
  ggplot(aes(x =phenotype, y = cuproptosis1, fill = phenotype)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.line.x = element_line(size = 0.8, linetype = "solid"),
        axis.line.y = element_line(size = 0.8, linetype = "solid"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20)  # 设置标题居中和字体大小
  ) +
  facet_wrap(~Celltype) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     size = 5,
                     label.y =1.2 ,
                     hide.ns = T) +
  scale_fill_manual(values = c("#33A02C","#FF7F00")) +
  labs(title = 'Allergen') +
  rotate_x_text(angle = 45) +
  ylab("Cuproptosis Score") +
  guides(fill = "none")
ggsave('./step1-scRNA/Allergen_cuproptosis.pdf',p6,height =8,width = 12,dpi = 500 )

p7 <- scRNA@meta.data[scRNA@meta.data$sample=='Dil',] %>%
  ggplot(aes(x =phenotype, y = cuproptosis1, fill = phenotype)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.line.x = element_line(size = 0.8, linetype = "solid"),
        axis.line.y = element_line(size = 0.8, linetype = "solid"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20)  # 设置标题居中和字体大小
  ) +
  facet_wrap(~Celltype) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     size = 5,
                     label.y =1.2 ,
                     hide.ns = T) +
  scale_fill_manual(values = c("#33A02C","#FF7F00")) +
  labs(title = 'Diluent') +
  rotate_x_text(angle = 45) +
  ylab("Cuproptosis Score") +
  guides(fill = "none")
ggsave('./step1-scRNA/Diluent_cuproptosis.pdf',p7,height =8,width = 12,dpi = 500 )

p8 <- scRNA@meta.data[scRNA@meta.data$sample=='Pre',] %>%
  ggplot(aes(x =phenotype, y = cuproptosis1, fill = phenotype)) +
  theme_pubr(base_size = 5) +
  theme(legend.position = "right", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.line.x = element_line(size = 0.8, linetype = "solid"),
        axis.line.y = element_line(size = 0.8, linetype = "solid"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20)  # 设置标题居中和字体大小
  ) +
  facet_wrap(~Celltype) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, lwd = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     size = 5,
                     label.y =1.2 ,
                     hide.ns = T) +
  scale_fill_manual(values = c("#33A02C","#FF7F00")) +
  labs(title = 'Baseline') +
  rotate_x_text(angle = 45) +
  ylab("Cuproptosis Score") +
  guides(fill = "none")
ggsave('./step1-scRNA/Baseline_cuproptosis.pdf',p8,height =8,width = 12,dpi = 500 )


genes_to_check <- cu$V1
p10 <- FeaturePlot(object = scRNA, features = genes_to_check, coord.fixed = T,pt.size = 0.1,
            cols = viridis(10),  order = T,ncol = 4)
ggsave('./step1-scRNA/FeaturePlot.pdf',p10,height =20,width = 16,dpi = 500 )
save(scRNA,file = './00.data/scRNA_Plot1.RData')

Idents(scRNA) <- 'Celltype'
scRNA1 <- subset(scRNA, idents = c("CD4 T cells"))
Idents(scRNA) <- 'sample'
scRNA1 <- subset(scRNA, idents = c("Ag"))
Idents(scRNA1) <- 'phenotype'
diff.mast <- FindAllMarkers(scRNA1)
diff.mast = diff.mast[diff.mast$cluster=='AA',]
gene <- cu$V1
diff1 = diff.mast %>% select(gene, everything()) %>% subset(p_val_adj<0.05&abs(avg_log2FC)>0.25)
diff_cu <- FindMarkers(scRNA1,ident.1 = 'AA',ident.2 = 'ANA',assay = 'RNA',slot = 'data')
cu$V1%in%rownames(diff_cu)


p100 <- FeaturePlot(scRNA1, features =gene, split.by = "phenotype", max.cutoff = 3,
            cols = c("grey", "red"))
ggsave('./step1-scRNA/FeaturePlot_all_AG.pdf',p100,height =48,width = 10,dpi = 500 )


FeaturePlot(scRNA, features = c("NLRP3", "GLS"), 
            blend = TRUE, 
            cols = c("gray80","red", "green"), 
            pt.size = 0.5, raster = F) +  
  theme(aspect.ratio = 1)
