library(openxlsx)
library(WGCNA)
library(ggplot2)
library(stringr)
library(pals)###我个人特别喜欢其中jet的配色方案
rm(list = ls())

gettreeplot<-function(hclustdata,groupdata,output,main){
  colordata<-data.frame(lapply(groupdata, function(i){
    c<-jet(length(unique(i))+2)[1:length(unique(i))+1]
    names(c)<-unique(i)
    return(sapply(i, function(j){c[j]}))
  }))
  rownames(colordata)<-rownames(groupdata)
  dhc <- as.dendrogram(hclustdata,hang=.01,frame.plot=F)
  colLab<-function(n){
    if(is.leaf(n)){
      a=attributes(n)
      labcol<-colordata[,1][rownames(colordata)==a$label]
      attr(n,"nodePar")<-c(a$nodePar,list(lab.col=labcol,lab.font=2,lab.cex=1,cex=0))
      if (ncol(colordata)>1){
        pointcol<-colordata[,2][rownames(colordata)==a$label]
        attr(n,"nodePar")<-c(a$nodePar,list(col=pointcol,cex=1.5))
      }
    }
    return(n)
  }
  
  lc<-unique(colordata[,1])
  if(ncol(colordata)>1){lc<-c(lc,unique(colordata[,2]))}
  ln<-unique(groupdata[,1])
  if(ncol(groupdata)>1){ln<-c(ln,unique(groupdata[,2]))}
  lp<-rep(4,length(unique(colordata[,1])))
  if(ncol(colordata)>1){lp<-c(lp,rep(20,length(unique(colordata[,2]))))}
  dL <- dendrapply(dhc, colLab)
  png(paste0(output,'.png'),width = 1000,height = 700)
  plot(dL,main=main)
  legend("topright", 
         legend = ln, col = lc,pch = lp, bty = "n",  pt.cex = 1.5, cex = 0.8 , 
         text.col = "black", horiz = FALSE, inset = c(0, 0.1))
  dev.off()
  pdf(paste0(output,'.pdf'))
  plot(dL , main=main)
  legend("topright", 
         legend = ln, col = lc,pch = lp, bty = "n",  pt.cex = 1.5, cex = 0.8 , 
         text.col = "black", horiz = FALSE, inset = c(0, 0.1))
  dev.off()
}



####多线程与内存设定
options(stringsAsFactors = F)
## 打开多线程
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=allowWGCNAThreads()-1
##设置内存上限为0.75*最大内存
memory.limit()
memory.limit(size = 0.75*(memory.limit()/(1024^3)))



### 数据读取与预处理
raw<-read.xlsx('C and without - with.xlsx',rowNames = T,check.names = F,sep.names = " ")
# rownames(raw)<-sapply(1:nrow(raw), function(i){paste0(raw[i,c('platform','id')],collapse = '_')})
data<-raw#[,grep('Control|Before|After',colnames(raw))]
# data=data[,-which(colnames(data)=='Control21')]
# rsd<-sapply(1:nrow(data), function(i){sd(data[i,])/mean(as.numeric(data[i,]))})
# data<-data[rsd>quantile(rsd,0.25),]
# vars<-sapply(1:nrow(data), function(i){var(as.numeric(data[i,]))})
cvs<-sapply(1:nrow(data), function(i){sd(as.numeric(data[i,]))/mean(as.numeric(data[i,]))})
data<-data[cvs>quantile(cvs,0.5),]
# annodata<-read.xlsx('../annotation.xlsx')
# rownames(annodata)<-sapply(1:nrow(annodata), function(i){paste0(annodata[i,c('platform','id')],collapse = '_')})
# annodata<-annodata[rownames(data),]


datExpr<-as.data.frame((t((data))))####WGCNA分析的输入数据，行为样本，列为基因
# datExpr=as.data.frame(scale(t((data))))####WGCNA分析的输入数据，行为样本，列为基因
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
## 表型形状文件
##行为样本，列为性状。性状以0,1代表有无
samples<-data.frame(type=unlist(sapply(colnames(data), function(i){str_split(i,'_')[[1]][1]})))
# samples$type<-sapply(samples$type,function(i){ifelse(i=='Without LN','Without LN',i)})
samples$type<-sapply(samples$type,function(i){ifelse(i=='Without LNC','Without LN',i)})
rownames(samples)<-colnames(data)
groupdata<-samples

for (ii in unique(samples$type)) {
  samples[,ii]<-as.numeric(samples$type==ii)
}
samples<-samples[,-1]
temps<-'WGCNA_Analysis'
if(!dir.exists(temps)){dir.create(temps)}

colordata<-data.frame(lapply(groupdata, function(i){
  c<-jet(length(unique(i))+2)[1:length(unique(i))+1]
  names(c)<-unique(i)
  return(sapply(i, function(j){c[j]}))
}))
rownames(colordata)<-rownames(groupdata)








#####样本聚类检查离群值##

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
temp0<-paste0(temps,'/1 Dendrogram')
if(!dir.exists(temp0)){dir.create(temp0)}
sampleTree = hclust(dist(scale(datExpr)), method = "average")
gettreeplot(sampleTree,groupdata,paste0(temp0,'/Dendrogram'),'Sample clustering to detect outliers')









####软阈值计算
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft<-pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
temp0<-paste0(temps,'/2 Soft_Threshold')
if(!dir.exists(temp0)){dir.create(temp0)}
write.csv(sft$fitIndices,paste0(temp0,'/Soft_Threshold.csv'),row.names = F)
png(paste0(temp0,"/Beta_Value.png"),width = 800,height = 600)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
pdf(paste0(temp0,"/Beta_Value.pdf"))
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



# ####一步直接构建加权共表达网络（Weight co-expression network)
# 
# ####构建加权共表达网络（Weight co-expression network)
# net = blockwiseModules(
#   datExpr,
#   power = sft$powerEstimate,
#   #maxBlockSize = 6000,
#   TOMType = "unsigned", minModuleSize = 30,
#   minClusterSize=30,#detectCutHeight = 0.995,
#   reassignThreshold = F, mergeCutHeight = 0.25,
#   numericLabels = TRUE, pamRespectsDendro = FALSE,
#   saveTOMs = F,deepSplit = 2,
#   verbose = 3#,indent = -2
# )
# table(net$colors) ###第一个是灰色模块基因，如果数量太多，则可能需要更改前面基因过滤步骤
# # Convert labels to colors for plotting
# # mergedColors1 = labels2colors(net$colors)
# # table(mergedColors1)
# moduleColors=mergedColors
# temp0<-paste0(temps,'/3 Weight_Co-expression_Network')
# if(!dir.exists(temp0)){dir.create(temp0)}
# # Plot the dendrogram and the module colors underneath
# png(paste0(temp0,"/Protein_Tree.png"),width = 800,height = 600)
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()
# ## assign all of the gene to their corresponding module 
# ## hclust for the genes.





####分步构建加权共表达网络（Weight co-expression network)
#(1)网络构建 Co-expression similarity and adjacency 
temp0<-paste0(temps,'/3 Weight_Co-expression_Network')
if(!dir.exists(temp0)){dir.create(temp0)}
adjacency = adjacency(datExpr, power = sft$powerEstimate) 
write.csv(adjacency,paste0(temp0,'/adjacency.csv'))
#(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
write.csv(TOM,paste0(temp0,'/TOM.csv'))
dissTOM = 1-TOM
# write.csv(dissTOM,paste0(temp0,'/dissTOM.csv'))
# (3) 聚类拓扑矩阵 Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
#(4) 聚类分支的修整 dynamicTreeCut 
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# 聚类树动态切割划分模块
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,cutHeight=0.997,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
#4. 绘画结果展示
# 模块名称由数字改为颜色
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(12,9)
png(paste0(temp0,"/Protein_Tree.png"),width = 800,height = 600)
# 聚类树
plot(geneTree, xlab="", sub="", main = "Protein clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
# 模块颜色
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colors")
dev.off()
pdf(paste0(temp0,"/Protein_Tree.pdf"))
# 聚类树
plot(geneTree, xlab="", sub="", main = "Protein clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
# 模块颜色
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colors")
dev.off()

#5. 聚类结果相似模块的融合，Merging of modules whose expression profiles are very similar
#在聚类树中每一leaf是一个短线，代表一个基因，
#不同分之间靠的越近表示有高的共表达基因，将共表达极其相似的modules进行融合
# 计算特征基因
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# 计算模块特征基因的相异性
MEDiss = 1-cor(MEs)
# 模块特征基因聚类
METree = hclust(as.dist(MEDiss), method = "average")
#选择有75%相关性的进行融合
MEDissThres = 0.25
# Plot the result
png(paste0(temp0,"/Module_Eigengenes_Tree.png"),width = 800,height = 600)
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
pdf(paste0(temp0,"/Module_Eigengenes_Tree.pdf"))
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
# 调用函数进行模块合并
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres)
# 模块颜色
mergedColors = merge$colors
# 新模块的特征基因:
mergedMEs = merge$newMEs








####模块与性状关联
temp0<-paste0(temps,'/4 Moudle_Train')
if(!dir.exists(temp0)){dir.create(temp0)}
moduleColors <- mergedColors
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) ##不同颜色的模块的ME值矩 (样本vs模块)
moduleTraitCor = cor(MEs, samples , use = "p")
write.xlsx(moduleTraitCor,paste0(temp0,'/corr.xlsx'),rowNames=T)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
write.xlsx(moduleTraitPvalue,paste0(temp0,'/Pvalue.xlsx'),rowNames=T)


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
png(paste0(temp0,"/Relationships.png"),width = 800,height = 600,res = 120)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(samples),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))
dev.off()
pdf(paste0(temp0,"/Relationships.pdf"))
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(samples),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))
dev.off()
###性状与模块显著性作图
if(!dir.exists(paste0(temp0,'/Trait'))){dir.create(paste0(temp0,'/Trait'))}
for (jj in colnames(samples)) {
  temp1<-paste0(temp0,'/Trait/',jj)
  if(!dir.exists(temp1)){dir.create(temp1)}
  trait<-as.data.frame(samples[,jj])
  names(trait)<-trait
  GeneSignificance=abs(as.numeric(cor(trait,datExpr, use="p")))
  #模块显著性为基因显著性的均值
  ModuleSignificance=tapply(GeneSignificance,
                            moduleColors, mean, na.rm=T)
  
  png(paste0(temp1,"/ Module-significance_with_",jj,'_bar.png'),width = 800,height = 600,res = 120)
  par(mfrow = c(1,1.5))
  plotModuleSignificance(GeneSignificance,moduleColors,ylim=c(0,1))
  dev.off()
  pdf(paste0(temp1,"/ Module-significance_with_",jj,'_bar.pdf'))
  plotModuleSignificance(GeneSignificance,moduleColors,ylim=c(0,1))
  dev.off()
  
  png(paste0(temp1,"/ Module-significance_with_",jj,'_box.png'),width = 800,height = 600,res = 120)
  par(mfrow = c(1,1.5))
  plotModuleSignificance(GeneSignificance,moduleColors,ylim=c(0,1),boxplot = T)
  dev.off()
  pdf(paste0(temp1,"/ Module-significance_with_",jj,'_box.pdf'))
  plotModuleSignificance(GeneSignificance,moduleColors,ylim=c(0,1),boxplot = T)
  dev.off()
  
}
##模块与性状关联作图
if(!dir.exists(paste0(temp0,'/Moudel'))){dir.create(paste0(temp0,'/Moudel'))}
mes<-as.matrix(mergedMEs)
write.xlsx(mes,paste0(temp0,'/Moudel/MEs.xlsx'),rowNames=T)
for (gg in colnames(mes)) {
  if (gg=='MEgrey'){next()}
  temp1<-paste0(temp0,'/Moudel/',gg)
  if(!dir.exists(temp1)){dir.create(temp1)}
  c<-unlist(sapply(rownames(mes), function(i){colordata[i,'type']}))
  p<-ggplot()+
    geom_bar(aes(y=mes[,gg],x=rownames(mes)),stat = 'identity',color=str_remove(gg,'ME'),fill=c)+
    xlab('Trait')+
    ylab('Moudle Eigengene Expression')+
    ggtitle(gg)+
    theme_bw()+
    theme(
      panel.background = element_rect(colour = 'white',fill ='white'),
      panel.grid = element_blank(),
      # panel.grid.major = element_blank(),
      # panel.border=element_blank(),
      legend.margin=margin(1.5,3,1.5,3,unit ="cm"),
      
      axis.ticks = element_line(size = 1.2,linetype=1),
      axis.line= element_line(size = 1.2,linetype=1),
      axis.ticks.length = unit(.31, "cm"),
      axis.text.x = element_text(face = "bold", color = c, size = 12,angle = 30,hjust = 1),
      axis.text.y = element_text(face = "bold", color = str_remove(gg,'ME'), size = 15),
      axis.title.x = element_text(face = "bold", color = "black", size = 20),
      axis.title.y = element_text(face = "bold", color = "black", size = 20),
      plot.title = element_text(face = "bold", color = str_remove(gg,'ME'), size = 28,hjust = 0.5),
      
      strip.text=element_text(face = "bold",angle=180,size=12),
      strip.background = element_rect(size=0.5,colour = 'white')
    )
  ggsave(paste0(temp1,'/',gg,' Eigengene-Expression.png'),p,height = 10,width = nrow(samples)*0.3,limitsize = F)
  ggsave(paste0(temp1,'/',gg,' Eigengene-Expression.pdf'),p,height = 10,width = nrow(samples)*0.3,limitsize = F)
  
}









###模块分析
###富集分析，GSEA分析的R脚本化 需要花点时间整理一下
###背景文件也要花时间整理一下
temp0<-paste0(temps,'/5 Moudle_Analysis')
if(!dir.exists(temp0)){dir.create(temp0)}
moudledata<-data.frame(mergedColors,row.names = colnames(datExpr))
# names (colors) of the modules
modNames = substring(names(MEs), 3)
#算出每个模块跟基因的皮尔森相关系数矩阵 (基因vs模块)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
names(geneModuleMembership) = paste("MM", modNames, sep="")
write.xlsx(cbind(moudledata,geneModuleMembership),paste0(temp0,'/moudle_cor.xlsx'),rowNames=T)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(MMPvalue) = paste("p.MM", modNames, sep="")
write.xlsx(cbind(moudledata,MMPvalue),paste0(temp0,'/moudle_pvalue.xlsx'),rowNames=T)


for (kk in unique(moudledata$mergedColors)) {
  if (kk=='grey'){next()}
  temp1<-paste0(temp0,'/',kk)
  if(!dir.exists(temp1)){dir.create(temp1)}
  ##模块内基因
  moudlegene<-rownames(moudledata)[moudledata$mergedColors==kk]
  
  ##模块内基因绘制聚类树
  ###发现对数据标量化，聚类效果反而会更好
  moudlesampleTree = hclust(dist(scale(datExpr[,moudlegene])), method = "average")
  gettreeplot(moudlesampleTree,groupdata,paste0(temp1,'/Dendrogram'),'Sample clustering to detect outliers')
  ##模块内基因绘制热图
  library(pheatmap)
  colanno<-colordata
  for (i in colnames(samples)) {
    colanno$type[samples[,i]==1]<-i
  }
  annocolor<-list(type=unique(colordata$type))
  names(annocolor$type)=unique(colanno$type)
  limits<-floor(min(max(t(scale(datExpr[,moudlegene]))),
                    -1*min(t(scale(datExpr[,moudlegene])))))
  bk <- c(seq(-1*limits,-0.1,by=0.01),seq(0,limits,by=0.01))
  c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
  pheatmap(t(scale(datExpr[,moudlegene])),show_rownames = F,annotation_colors = annocolor,fontsize_col=5,
              annotation_col=colanno,border_color = NA,color = c,breaks=bk,angle_col=90,
              filename = paste0(temp1,'/',kk,'.pdf'))
  pheatmap(t(scale(datExpr[,moudlegene])),show_rownames = F,annotation_colors = annocolor,
           annotation_col=colanno,border_color = NA,color = c,breaks=bk,angle_col=90,fontsize_col=5,
           filename = paste0(temp1,'/',kk,'.png'))
  
  ##GS-MM散点图
  for (jj in colnames(samples)) {
    temp2<-paste0(temp1,'/',jj)
    if(!dir.exists(temp2)){dir.create(temp2)}
    
    trait<-as.data.frame(samples[,jj])
    names(trait)<-trait
    geneTraitSignificance=as.data.frame(cor(datExpr, trait, use = "p"))
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
    
    names(geneTraitSignificance) = paste("GS.", jj, sep="")
    names(GSPvalue) = paste("p.GS.", jj, sep="")
    gsmm<-data.frame(kk=abs(geneModuleMembership[moduleColors==kk, paste0('MM',kk)]),
                     jj=abs(geneTraitSignificance[moduleColors==kk, 1]))
    colnames(gsmm)<-c(kk,jj)
    rownames(gsmm)<-rownames(geneTraitSignificance)[moduleColors==kk]
    write.xlsx(gsmm,paste0(temp2,'/',kk,'.xlsx'),rowNames=T)

    png(paste0(temp2,'/GS-MM.png'),width = 800,height = 600)
    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleColors==kk, paste0('MM',kk)]),
                       abs(geneTraitSignificance[moduleColors==kk, 1]),
                       xlab = paste("Module Membership in", kk, "module"),
                       ylab = "Gene significance for Luminal",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = kk)
    dev.off()
    pdf(paste0(temp2,'/GS-MM.pdf'))
    #sizeGrWindow(7, 7);
    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleColors==kk, paste0('MM',kk)]),
                       abs(geneTraitSignificance[moduleColors==kk, 1]),
                       xlab = paste("Module Membership in", kk, "module"),
                       ylab = "Gene significance for Luminal",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = kk,fill=kk)
    dev.off()
  }
 
  
}








####模块内hubgene分析
#HubGenes，连接度top10%  模块100-1000个基因
#模块基因<100,top10
#模块基因>1000，top100
temp0<-paste0(temps,'/6 HubGene_Analysis')
if(!dir.exists(temp0)){dir.create(temp0)}
# hubgenes<-data.frame(matrix(nrow=nrow(datExpr),ncol = 1))
# rownames(hubgenes)<-rownames(datExpr)
hubgenes<-c()
for (kk in unique(moudledata$mergedColors)) {
  if (kk=='grey'){next()}
  temp1<-paste0(temp0,'/',kk)
  if(!dir.exists(temp1)){dir.create(temp1)}
  ##模块内基因
  moudlegene<-rownames(moudledata)[moudledata$mergedColors==kk]
  ##top 10~1000
  top<-0.1*length(moudlegene)
  if (top<10) {
    top<-10
  }else if(top>100){
    top<-100
  }
  ###连接度
  IMConn = softConnectivity(datExpr[, moudlegene])
  write.xlsx(data.frame(softConnectivity=IMConn,row.names = moudlegene),
             paste0(temp1,'/',kk,'_softConnectivity.xlsx'),rowNames=T)
  hub<-(rank(-IMConn) <= top)
  hubgene<-datExpr[,moudlegene][hub]
  write.xlsx(hubgene,paste0(temp1,'/',kk,'_hubgene.xlsx'),rowNames=T)
  # hubgenes<-cbind(hubgenes,hubgene)
  hubgenes<-c(hubgenes,colnames(hubgene))
  names(hubgenes)[length(hubgenes):(length(hubgenes)+ncol(hubgene))-ncol(hubgene)]<-rep(kk,ncol(hubgene))
  ###hubgene绘制热图
  limits<-floor(min(max(t(scale(hubgene))),-1*min(t(scale(hubgene)))))
  bk <- c(seq(-1*limits,-0.1,by=0.01),seq(0,limits,by=0.01))
  c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
  pheatmap(t(scale(hubgene)),show_rownames = F,annotation_colors = annocolor,fontsize_col=5,
           annotation_col=colanno,border_color = NA,color = c,breaks=bk,angle_col=90,
           filename = paste0(temp1,'/',kk,'.pdf'))
  pheatmap(t(scale(hubgene)),show_rownames = F,annotation_colors = annocolor,
           annotation_col=colanno,border_color = NA,color = c,breaks=bk,angle_col=90,fontsize_col=5,
           filename = paste0(temp1,'/',kk,'.png'))
}
select<-colnames(datExpr) %in% hubgenes
selectTOM = dissTOM[select, select]
plotDiss = selectTOM^20
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
png( paste0(temp0,'/HubGene_dissTOM.png'),width = 800,height = 600)
TOMplot(plotDiss, selectTree, selectColors, main = "dissTOM^20 heatmap plot, hungene genes")
dev.off()
pdf( paste0(temp0,'/HubGene_dissTOM.pdf'))
TOMplot(plotDiss, selectTree, selectColors, main = "dissTOM^20 heatmap plot, hungene genes")
dev.off()

##所有HunGene绘制热图
#hubgenes<-hubgenes[,-1]
hubgenes<-datExpr[,hubgenes]
limits<-floor(min(max(t(scale(hubgenes))),-1*min(t(scale(hubgenes)))))
bk <- c(seq(-1*limits,-0.1,by=0.01),seq(0,limits,by=0.01))
c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
pheatmap(t(scale(hubgenes)),show_rownames = F,annotation_colors = annocolor,fontsize_col=5,
         annotation_col=colanno,border_color = NA,color = c,breaks=bk,angle_col=90,
         filename = paste0(temp0,'/HubGene.png'))
pheatmap(t(scale(hubgenes)),show_rownames = F,annotation_colors = annocolor,
         annotation_col=colanno,border_color = NA,color = c,breaks=bk,angle_col=90,fontsize_col=5,
         filename = paste0(temp0,'/HubGene.pdf'))





###绘制网络图
###是否考虑可以和hubgene合并呢！！
###show自定义取值，all为所有模块基因
shows='all'
library(flashClust)
library(igraph)
temp0<-paste0(temps,'/9 Network_Analysis')
if(!dir.exists(temp0)){dir.create(temp0)}
for (kk in unique(moudledata$mergedColors)) {
  if (kk=='grey'){next()}
  temp1<-paste0(temp0,'/',kk)
  if(!dir.exists(temp1)){dir.create(temp1)}
  inModule<-moudledata$mergedColors==kk
  modTOM<-TOM[inModule, inModule]
  if (shows=='all') {show=sum(inModule)}
  moudlegene<-rownames(moudledata)[moudledata$mergedColors==kk]
  ###连接度
  IMConn = softConnectivity(datExpr[, moudlegene])
  show<-(rank(-IMConn) <= show)
  cyt = exportNetworkToCytoscape( modTOM[show,show],
                                  edgeFile = paste(temp1,"/CytoscapeInput-edges-", kk, ".txt", sep=""),
                                  nodeFile = paste(temp1,"/CytoscapeInput-nodes-", kk, ".txt", sep=""),
                                  weighted = TRUE,
                                  threshold = 0.05,
                                  nodeNames = moudlegene[show],
                                  altNodeNames = moudlegene[show],
                                  nodeAttr = moduleColors[inModule][show])
  ######################################绘制基因网络草图############################
  Edge = read.table(paste(temp1,"/CytoscapeInput-edges-", kk, ".txt", sep=""),header = 1)
  # Etop = (rank(-Edge['weight'])<=200)
  # Edge = Edge[Etop,]
  gd = graph.data.frame(Edge[,c(1,2,1)],directed = F)

  for (ee in E(gd)) {
    
  }
  pdf(paste0(temp1,'/NetGraph.pdf'),width = 15,height = 15,onefile=FALSE);
  plot(gd,
       layout=layout.fruchterman.reingold
)
  dev.off()
}
