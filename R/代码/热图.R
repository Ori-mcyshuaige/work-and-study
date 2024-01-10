##如果没有R包先下载R包
# install.packages("openxlsx")
# install.packages("pheatmap")
# install.packages("pals")

##导入R包
library(openxlsx)
library(pheatmap)
library(pals)
##建立数据所在文件路径
setwd("D:/马驰宇/桌面/何院专辑/研究院项目/戴主任项目重分析/氧化脂质/新数据补值/氧化脂质/加入临床指标/热图")
##导入所需数据
data<-read.xlsx("785个蛋白.xlsx",sheet = 1,colNames = T,rowNames = T,check.names = F,sep.names = " ")
# annotation_col<-read.xlsx("氧化脂质与抗体表达之间的关系.xlsx",sheet = 3,colNames = T,rowNames = T,check.names = F,sep.names = " ")
# annotation_row<-read.xlsx("显著性检验Control-active-inactive.xlsx",sheet = 4,colNames = T,rowNames = T,check.names = F,sep.names = " ")
# display_number<-read.xlsx("显著性检验Control-active-inactive.xlsx",sheet = 4,colNames = T,rowNames = T,check.names = F,sep.names = " ")
# data<-data[,rownames(annotation_col)]
# colnames(data)<-rep(c("CONT","Model01","SM"),each=6)


df <- read.csv('./01.差异分析/df_diff_GSE98793.csv',header = T,row.names = 1)
# data <- read.csv('./01.差异分析/ori_GSE98793.csv',header = T,row.names = 1)
# data <- read.csv('./01.差异分析/batch_GSE98793.csv',header = T,row.names = 1)
annotation_col <- as.data.frame(str_split_i(colnames(data),'_',1))
rownames(annotation_col) <- colnames(data)
colnames(annotation_col) <- 'Group'
data <- df[1:104,]
for(i in colnames(data)){
  data[,i] <- as.numeric(data[,i])
}
#数据转换
#对行做z-score
for(ii in 1:nrow(data)){
      result<-sapply(as.numeric(data[ii,]),function(i){(i-mean(as.numeric(data[ii,])))/sd(as.numeric(data[ii,]))})
      data[ii,]<-result
}
#log转化，缩小数据间距，调整色阶，原则上先不做
for(ii in rownames(data)){
  result<-sapply(as.numeric(data[ii,]),function(i){log10(i)})
  data[ii,]<-result
}

##注释配色
sample<-jet(length(unique(annotation_col[,1]))+3)[1:length(unique(annotation_col[,1]))+1]##jet()函数可以用于挑选颜色，例：jet(5)可以挑选5个颜色
names(sample)<-unique(annotation_col[,1])##给要添加的注释分配颜色
meta1<-jet(length(unique(annotation_row[,1])))
names(meta1)<-unique(annotation_row[,1])
# meta2<-jet(length(unique(annotation_row[,2])))
# names(meta2)<-unique(annotation_row[,2])
annotation_color<-list(group=sample,"Classification I"=meta1)##将添加的注释创建一个列表储存起来，作为pheatmap的参数
data<-data[,c(grep("Control",annotation_col[,1]),grep("Inactive",annotation_col[,1]),c(1:227)[-c(grep("Inactive",annotation_col[,1]),grep("Control",annotation_col[,1]))])]
data <- data[,c(grep("control",annotation_col[,1]),grep("MDD",annotation_col[,1]))]
##热图每个格子的颜色
limits<-floor(max(max(data,na.rm = T),1*min(data,na.rm = T))*100000)/100000
bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))


##pheatmap函数绘制热图
pheatmap(data,scale = "none",color = c,show_colnames=T,show_rownames = F,cluster_cols=F,cluster_rows = T,annotation_col = annotation_col,#clustering_method = "complete",#annotation_row = annotation_row,
         breaks=1*bk,cellwidth = 20, cellheight = 1, fontsize_row  = 10, fontsize_col = 10, border_color = "white",annotation_colors = annotation_color,
         filename = "heatmap.pdf")
pheatmap(data,scale = "none",color = c,show_colnames=F,cluster_cols=F,cluster_rows = F,#annotation_col = annotation_col,annotation_row = annotation_row,
         breaks=1*bk,cellwidth = 30, cellheight = 10, fontsize_row  = 9, fontsize_col = 3, border=FALSE,#annotation_colors = annotation_color,
         filename = "heatmap.pdf")


pheatmap(data,color=c,breaks=0.24*bk,annotation_col = annotation_col,annotation_row = annotation_row,border=F,
         display_numbers = display_number,fontsize_number = 2,cluster_cols=F,cluster_rows = F,show_colnames=F,show_rownames=F,
         annotation_colors = annotation_color,cellwidth = 4, cellheight = 7, fontsize_row  = 4, fontsize_col = 4,
         filename = "heatmap.png")

