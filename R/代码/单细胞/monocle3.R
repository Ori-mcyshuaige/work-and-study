# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
# dir.create("Monocle3")



# #导入数据
# #在这里你需要导入三个数据：细胞的表达矩阵，meta信息，gene_annotation
# data <- 细胞的表达矩阵
# cell_metadata <- meta信息
# gene_annotation <- data.frame(gene_short_name = rownames(data))
# rownames(gene_annotation) <- rownames(data)
# 
# #创建CDS对象
# cds <- new_cell_data_set(data,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)  #这样就获得了monocle的对象了，可类比于创建seurat对象
# 
# #思考一下我们在创建了seurat对象后要干嘛呢？降维聚类分群，那么monocle也不例外
# 
# #NormalizeData+ScaleData+RunPCA
# cds <- preprocess_cds(cds, num_dim = 50)     #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
# plot_pc_variance_explained(cds)    #像seurat一样展示pc数
# 
# #umap降维
# cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
# 
# #tSNE降维
# cds <- reduce_dimension(cds, reduction_method="tSNE")
# 
# #聚类
# cds <- cluster_cells(cds)   
# 
# #那么我们到这里就完成了分群了，后面就是对于每一个群的细胞注释和细胞类型的结果展示了，这个大家可以去看周运来老师的简书，而我的目的不是这个
scRNA <- readRDS('./final_scRNA_rpca.rds')

##创建CDS对象并预处理数据
#首先我们要得到seurat对象的表达矩阵
data <- GetAssayData(scRNA, assay = 'RNA', slot = 'counts')

#其次获得meta信息
cell_metadata <- scRNA@meta.data

#然后要获得gene_annotation 
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

#创建CDS对象
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#然后就是降维聚类分群
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)     #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
plot_pc_variance_explained(cds)    #像seurat一样展示pc数

#umap降维
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA

#tSNE降维
cds <- reduce_dimension(cds, reduction_method="tSNE")

#聚类
cds <- cluster_cells(cds,cluster_method = 'louvain') 

##是否需要把monocle分析的umap结果替换为seurat的umap结果？？？？
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed   #因此我们以后只要是画umap的点图就跟seurat的图片形状是一样的

#画图看一下
plot_cells(cds, reduction_method="UMAP", color_cells_by="Cell.types") 

#######拟时序分析

#第一步:轨迹推断
cds <- learn_graph(cds)     #是的，就这么简单，用 learn_graph即可

#画图  我们可以根据meta信息来画出相关的图片
head(colData(cds))
plot_cells(cds,
           color_cells_by = "Cell.types",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1.5)                                    #这个图就可以看出细胞直接的分化轨迹了
#黑色的线显示的是graph的结构。
#数字带白色圆圈表示不同的结局，也就是叶子。
#数字带黑色圆圈代表分叉点，从这个点开始，细胞可以有多个结局。
#这些数字可以通过label_leaves和label_branch_points参数设置。


##第二步：选择根 就是选择细胞的起源 这也是monocle3与monocle2最大的区别，手动选择(自定义)根 ，使用其中一种就可以！！！
#第一种：
cds = order_cells(cds)  #在前阵子有bug 无法展示，现在好像又可以了
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)

#第二种：用辅助线
p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(0,1,0.25))
embed <- data.frame(Embeddings(scRNA, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -6.75 & UMAP_1 < -6.5 & UMAP_2 > 0.24 & UMAP_2 < 0.25)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)    #这里就选择了根
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)

#第三种：定义一个函数
# a helper function to identify the root principal points:
get_earliest_principal_node  <- function(cds, time_bin="epithelial cells"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))    #很多人这里一直在哭诉error，那是都没有理解这一步在干嘛，很多无脑运行别人的代码，别人代码选择XX细胞作为根，你的数据集里又没有，所以报错说没有node啦

###差异分析
#计算基因按照轨迹的变化
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

#基因表达趋势图
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
#在这里运用到文章的话，应该可以这样用，你可以用这些变化top基因来讲一些东西，或者你可以show你自己fouce on的gene


##然后就是还可以做一个寻找共表达模块，怎么应用呢，就可以说是这个模块导致他们状态发生变化，然后查看一个模块里的基因，再讲一些东西

##寻找共表达模块
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$celltype)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")


####最后，我们可以把monocle的拟时序结果导入seurat对象里了
#提取拟时分析结果返回seurat对象
pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(scRNA@meta.data)]
scRNA$pseudotime <- pseudotime
p = FeaturePlot(scRNA, reduction = "umap", features = "pseudotime")
# pseudotime中有无限值，无法绘图。
ggsave("Pseudotime_Seurat.pdf", plot = p, width = 8, height = 6)
saveRDS(scRNA, file = "sco_pseudotime.rds")
