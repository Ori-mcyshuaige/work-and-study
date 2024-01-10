library(openxlsx)
library(pheatmap)
library(pals)

data <- as.data.frame(t(read.csv('05.通路富集分析/mean_dfrhygene.csv',row.names = 1,check.names = F)))

annotation_col <- as.data.frame(str_split_i(colnames(data),'_',1))
rownames(annotation_col) <- colnames(data)
colnames(annotation_col) <- 'Group'

amp <- read.csv('02.昼夜振幅估计/allampstat.csv',row.names = 1)

annotation_row <- data.frame(Amplitude=amp[rownames(data),4])
rownames(annotation_row) <- rownames(data)
annotation_row[,1] <- ifelse(annotation_row[,1]<0,'Higher in Ctrl','Higher in DM')
for(ii in 1:nrow(data)){
  result<-sapply(as.numeric(data[ii,]),function(i){(i-mean(as.numeric(data[ii,])))/sd(as.numeric(data[ii,]))})
  data[ii,]<-result
}

sample<-c('black','red')
names(sample)<-unique(annotation_col[,1])
gene<-jet(length(unique(annotation_row[,1]))+3)[1:length(unique(annotation_row[,1]))+1]##jet()函数可以用于挑选颜色，例：jet(5)可以挑选5个颜色
names(gene)<-unique(annotation_row[,1])
annotation_color <- list(Group=sample,Amplitude=gene)


limits<-floor(max(max(data,na.rm = T),1*min(data,na.rm = T))*100000)/100000
bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

data1 <- data[rownames(annotation_row)[annotation_row$Amplitude=='Higher in Ctrl'],
             c('Ctrl_ZT0','Ctrl_ZT4','Ctrl_ZT8','Ctrl_ZT12','Ctrl_ZT16','Ctrl_ZT20','DM_ZT0','DM_ZT4','DM_ZT8','DM_ZT12','DM_ZT16','DM_ZT20')]
data2 <- data[rownames(annotation_row)[annotation_row$Amplitude!='Higher in Ctrl'],
              c('Ctrl_ZT0','Ctrl_ZT4','Ctrl_ZT8','Ctrl_ZT12','Ctrl_ZT16','Ctrl_ZT20','DM_ZT0','DM_ZT4','DM_ZT8','DM_ZT12','DM_ZT16','DM_ZT20')]
data1 <- data1[order(sapply(1:nrow(data1),function(x){mean(unlist(data1[x,1:3]))})),]
data2 <- data2[order(sapply(1:nrow(data2),function(x){mean(unlist(data2[x,1:3]))})),]
dat <- rbind(data1,data2)
##pheatmap函数绘制热图
pheatmap(dat,scale = "none",color = c,show_colnames=T,show_rownames = F,cluster_cols=F,cluster_rows = F,annotation_col = annotation_col,annotation_row = annotation_row,#clustering_method = "complete",
         breaks=1*bk,cellwidth = 30, cellheight = 1, fontsize_row  = 10, fontsize_col = 10, border_color = "white",annotation_colors = annotation_color,
         filename = "heatmap.pdf")
pheatmap(data1,scale = "none",color = c,show_colnames=T,show_rownames = F,cluster_cols=F,cluster_rows = T,annotation_col = annotation_col,annotation_row = annotation_row,#clustering_method = "complete",
         breaks=1*bk,cellwidth = 30, cellheight = 1, fontsize_row  = 10, fontsize_col = 10, border_color = "white",annotation_colors = annotation_color,
         filename = "Higher in Ctrl heatmap.pdf")
pheatmap(data2,scale = "none",color = c,show_colnames=T,show_rownames = F,cluster_cols=F,cluster_rows = T,annotation_col = annotation_col,annotation_row = annotation_row,#clustering_method = "complete",
         breaks=1*bk,cellwidth = 30, cellheight = 1, fontsize_row  = 10, fontsize_col = 10, border_color = "white",annotation_colors = annotation_color,
         filename = "Higher in DM heatmap.pdf")
