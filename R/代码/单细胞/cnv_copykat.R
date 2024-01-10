library(Seurat)
library(copykat)
library(tidyverse)
library(infercnv)

# library(devtools)
# install_github("navinlabcode/copykat")

# library(rjags)
# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
# ####linux中安装先运行如下代码
# #####sudo apt-get install jags
# 
# BiocManager::install("rjags")
# BiocManager::install("infercnv",force = TRUE)

#####  copykat   ######

dir.create('04.pre_cnv')

counts <- as.matrix(data.combined_all@assays$RNA@counts)
ckat_rs <- copykat(rawmat = counts,ngene.chr = 5,sam.name = 'hcc_cir',n.cores = 30,
                   genome = "hg20",distance = "euclidean")
saveRDS(ckat_rs, "./04.pre_cnv/ckat_rs.rds")

mallignant <- read.delim("hcc_cir_copykat_prediction.txt", row.names = 1)
# mallignant <- data.frame('copykat.pred' = mallignant[rownames(data.combined_all@meta.data),],
#                          row.names = rownames(data.combined_all@meta.data))
data.combined_all <- AddMetaData(data.combined_all, metadata = mallignant)

p1 <- DimPlot(data.combined_all, group.by = "Cell.types", label = T,raster=FALSE)
p2 <- DimPlot(data.combined_all, group.by = "copykat.pred",raster=FALSE) + scale_color_manual(values = c("red", "gray",'black'))
p3 <- DimPlot(data.combined_all, group.by = "group", label = T,raster=FALSE)
pc <- p1 + p2 + p3
ggsave("./04.pre_cnv/pred_copykat.pdf", pc, width = 16, height = 5)



#####  infercnv   ####
inferscRNA <- subset(data.combined_all,subset = Cell.types=='Hepatocytes')
##第一个文件:单细胞RNA-seq表达量的原始矩阵
dfcount = as.data.frame(inferscRNA@assays$RNA@counts)
##第二个文件:注释文件，记录肿瘤和正常细胞
groupinfo= data.frame(cellId = colnames(dfcount),cellType= inferscRNA@meta.data$group )
##第三个文件:基因注释文件
library(AnnoProbe)
geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 

write.table(dfcount ,file ='./04.pre_cnv/expFile.txt',sep = '\t',quote = F)
write.table(groupinfo,file = './04.pre_cnv/metaFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file ='./04.pre_cnv/geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./04.pre_cnv/expFile.txt",
                                    annotations_file="./04.pre_cnv/metaFiles.txt",
                                    delim="\t",
                                    gene_order_file= "./04.pre_cnv/geneFile.txt",
                                    ref_group_names='cirrhotic') # 如果有正常细胞的话，把正常细胞的分组填进去
future::plan("multisession",workers=30)# 多核并行处理
options(future.globals.maxSize = 35 * 1000 * 1024^2)
options(scipen = 100)###Error in m[match(oldnodes, m)] <- 1:(N - 1) : NAs are not allowed in subscripted assignments
options(bitmapType="Xlib")###Error in plot.new() : figure margins too large
options(bitmapType='cairo')### x11
###ulimit -s unlimited  #c pan limit
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='./04.pre_cnv', 
                             cluster_by_groups=FALSE,  # 是否根据细胞注释文件的分组对肿瘤细胞进行
                             cluster_references = TRUE,
                             denoise=TRUE, #去噪
                             HMM=FALSE, #是否基于HMM预测CNV,选择F会加快运算速度
                             output_format = "pdf",
                             useRaster = FALSE)

infercnv.dend <- read.dendrogram(file = './04.pre_cnv/infercnv_subclusters.observations_dendrogram.txt')