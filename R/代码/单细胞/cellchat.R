# devtools::install_github("sqjin/CellChat")
library(Seurat)
##install.packages("remotes")
## remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(future)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 1000*1024^4)


scRNA <- readRDS('./final_scRNA_rpca.rds')

##创建cellchat对象
cellchat <- createCellChat(scRNA@assays$RNA@data,
                          meta = scRNA@meta.data,
                          group.by = 'Cell.types')
##查看每个cluster有多少个细胞
groupSize <- as.numeric(table(cellchat@idents))
groupSize

##导入配体受体数据库，有人有鼠
CellChatDB <- CellChatDB.human
str(CellChatDB)##查看数据库信息
##包含interaction/complex/cofactor/geneInfo四个dataframe
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() #下一步不出图的时候运行
showDatabaseCategory(CellChatDB)


unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object

###对表达数据进行预处理。首先在一个细胞组中识别过表达的配体或受体，
###然后将基因表达数据投射到PPI网络上，如果配体和受体过表达则识别配体和受体之间的
###相互作用。

## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)##multisession   ##multicore
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat <- identifyOverExpressedInteractions(cellchat)
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project

####推断细胞通讯网络
#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "./06.cellchat/net_lr.csv")

###推断信号通路水平的细胞通讯网络（结果存储在@netP下，有一个概率值和对应的pval）
###我们可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络。
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "./06.cellchat/net_pathway.csv")

##细胞互作关系展示
#所有细胞群总体观：细胞互作数量与强度统计分析
#统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,  weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")
# save as TIL_net_number_individual.pdf
##左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。发出箭头的细胞表达配体，箭头指向的细胞表达受体。配体-受体对越多，线越粗。右图：互作的概率/强度值（强度就是概率值相加）


##检查每种细胞发出的信号
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,  vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as ctype_TIL_net_number_individual.pdf

## 运行上述代码出现报错 Error in plot.new() : figure margins too large
# par("mar")
## [1] 5.1 4.1 4.1 2.1
# par(mar=c(1,1,1,1))
# 重新运行上面的代码
##每个细胞如何跟别的细胞互作（number of interaction图）

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as ctype_TIL_net_strength_individual.pdf
#每个细胞如何跟别的细胞互作（互作的强度/概率图）

####单个信号通路或配体-受体介导的细胞互作可视化（层次图、网络图、和弦图、热图）
cellchat@netP$pathways#查看都有哪些信号通路
# 选择其中一个信号通路，比如说TGFb
pathways.show <- c("TGFb")  
##层次图（Hierarchy plot）
levels(cellchat@idents)# show all celltype
vertex.receiver = c(1,2,4,6)## define a numeric vector （淋系细胞）giving the index of the celltype as targets
#par(mar=c(5.1,4.1,4.1,2.1))
netVisual_aggregate(cellchat, signaling = pathways.show,vertex.receiver = vertex.receiver)
# save as TIL/CXCL_hierarchy.pdf

###网络图（Circle plot）
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# save as TIL_CXCL_circle.pdf

#和弦图（Chord diagram）
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# save as TIL_CXCL_chord.pdf

#热图（Heatmap）
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# save as TIL_CXCL_heatmap.pdf

##配体-受体层级的可视化（计算各个ligand-receptor pair对信号通路的贡献）
#计算配体受体对选定信号通路的贡献值（在这里就是查看哪条信号通路对TGFb贡献最大）
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)#提取对TGFb有贡献的所有配体受体 
# save as TIL_CXCL_LR_contribution.pdf
##层次图（ Hierarchy plot）
#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.TGFb[1,] 
vertex.receiver = c(1,2,4,6) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# save as TIL_CXCL_hierarchy2.pdf

#网络图（Circle plot）
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# save as TIL_CXCL_circle2.pdf

#和弦图（Chord diagram）
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# save as TIL_CXCL_chord2.pdf


##自动（批量）保存每个信号通路的互作结果
# Access all the signaling pathways showing significant communications将所有信号通路找出来
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1,2,4,6) #不画层次图就不需要这一步
dir.create("./06.cellchat/all_pathways_com_circle")#创建文件夹保存批量画图结果
setwd("./06.cellchat/all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
            vertex.receiver = vertex.receiver, layout = "circle")#绘制网络图
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i],"_L-R_contribution.pdf"), 
         plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}

###多个配体-受体介导的细胞互作关系可视化
##气泡图（全部配体受体）
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9),
                     targets.use = c(1,2,4,6), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble.pdf", p, width = 8, height = 12)#髓系对淋巴的调节
# save as TIL/Mye_Lymph_bubble.pdf

##气泡图（指定信号通路或配体-受体）
#比如制定CCL和CXCL这两个信号通路
netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), targets.use = c(1,2,4,6), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)
##气泡图（指定信号通路或配体-受体并指定细胞）
# show all the significant interactions (L-R pairs) based on user's input
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","TGFb"))
netVisual_bubble(cellchat, sources.use = c(3,6,8), targets.use = c(1,4,5), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)

##参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）

## Plot the signaling gene expression distribution
p = plotGeneExpression(cellchat, signaling = "TGFb",raster=FALSE)
ggsave("TGFb_GeneExpression_vln.pdf", p, width = 8, height = 8)
p = plotGeneExpression(cellchat, signaling = "TGFb", type = "dot",color.use = 'blue')
ggsave("TGFb_GeneExpression_dot.pdf", p, width = 8, height= 6)

