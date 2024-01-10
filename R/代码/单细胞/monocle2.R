options(stringsAsFactors = F)
gc()
library(monocle)
library(Seurat)
library(dplyr)
library(patchwork)
library(celldex)
library(CytoTRACE)

scRNA <- readRDS('./final_scRNA_rpca.rds')

###挑选出我们感兴趣的细胞类型进行下一步分析，值得注意的是，
###做单细胞拟时序分析只能选择一类细胞或者是一类细胞的不同亚群

sc_BP <- scRNA[,scRNA$Cell.types %in% c('Hepatocytes')]####我们挑选了B细胞和浆细胞作为研究对象

###上述操作形成一个新的seurat对象，要对该对象重新进行降维、分群
sc_BP <- FindVariableFeatures(sc_BP,selection.method = 'vst',nfeatures = 2000)
sc_BP <- ScaleData(sc_BP)
sc_BP <- RunPCA(sc_BP,npcs = 30,verbose = F)
sc_BP <- RunUMAP(sc_BP,reduction = 'pca',dims = 1:30)
Idents(sc_BP) <- 'orig.ident'
DimPlot(sc_BP,reduction = 'umap')



####分群并绘图
sc_BP <- FindNeighbors(sc_BP,reduction = 'pca',dims = 1:30)
sc_BP <- FindClusters(sc_BP,resolution = 0.5)
DimPlot(sc_BP,group.by = 'RNA_snn_res.0.5')

####使用monocle2对上述seurat对象进行分析
expr_matrix <- as(as.matrix(sc_BP@assays$RNA@counts),'sparseMatrix')###提取稀疏矩阵
p_data <- sc_BP@meta.data
p_data$Cell.types <- sc_BP@active.ident###整合每个细胞的下拨鉴定信息到p_data
f_data <- data.frame(gene_short_name = row.names(sc_BP),row.names = row.names(sc_BP))###提取基因名称

pd <- new('AnnotatedDataFrame',data = p_data)
fd <- new('AnnotatedDataFrame',data = f_data)
####AnnotatedDataFrame对象包含测量变量及其元数据描述的类；元数据指的是对象的相关数据，关于对象的一切相关信息

######创建monocle2需要的cds对象
cds <- newCellDataSet(expr_matrix,#稀疏矩阵
                      phenoData = pd,#元数据，seurat对象中的meta.data
                      featureData = fd,#是基因名称
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)###标准化细胞细胞之间的mRNA差异
cds <- estimateDispersions(cds)###离散度值帮助我们进行下一步分析

cds <- detectGenes(cds = cds,min_expr = 0.1)##这一步可以在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))

#####过滤低质量细胞
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >=10))
#####step1.筛选差异表达基因，用于轨迹定义基因
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = '~Cell.types',cores = 1)#~后面表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)
deg <- subset(diff,qval<0.01)
deg <- deg[order(deg$qval,decreasing = F),] #排序

ordergene <- rownames(deg)
cds <- setOrderingFilter(cds,ordergene)#该函数标记将在随后对clusterCells的调用中用于聚类的基因

####step2.降维
cds <- reduceDimension(cds,max_components = 2,reduction_method = "DDRTree")


###step3.拟时间轴轨迹构建和在拟时间内排序细胞
cds1 <- orderCells(cds)


####step4.monocle2可视化拟时序
plot_cell_trajectory(cds1,color_by = 'Pseudotime',size=1,show_backbone = T)
plot_cell_trajectory(cds1,color_by = 'State',size=1,show_backbone = T)


####单单使用monocel2无法半段判断细胞发育起点，需要借助CytoTRACE包进一步分析
####利用monocle2创建的cds对象，用CytoTRACE包做拟时序分析
exp1 <- as.matrix(sc_BP@assays$RNA@counts)#提取稀疏矩阵
exp1 <- exp1[apply(exp1>0,1,sum)>5,] #过滤
####计算每个细胞的而分化程度数据
results <- CytoTRACE(exp1,ncores = 1)

mono_meta <- data.frame(t(cds1@reducedDimS),#####mono结果体现在这里
                        cds1$Pseudotime,
                        cds1$State,
                        cds1$Cell.types)
###获取元数据
head(mono_meta)

colnames(mono_meta) <- c('C1','C2','Pseudotime','State','Cell.types')
phenot <- as.character(cds$Cell.types)
names(phenot) <- rownames(mono_meta)
emb <- mono_meta[,1:2]

####可视化
plotCytoTRACE(results,phenoCell.types = phenot,emb = emb)
###对比两个包的图发现：CytoTRACE以左下角为发育起点，而monocel2以右下角为发育起点，
###在这里我们应该以CytoTRACE为准，因此需要更改monocle2的发育起点

cds2 <- orderCells(cds1,root_state = 3)

plot_cell_trajectory(cds2,color_by = 'Pseudotime',size=1,show_backbone = T)