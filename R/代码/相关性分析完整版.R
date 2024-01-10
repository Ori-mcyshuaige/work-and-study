library(openxlsx)
library(pheatmap)
library(fdrtool)
library(ggplot2)
library(networkD3)
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)
library(pals)

# setwd("D:/马驰宇/桌面/何院专辑/研究院项目/戴主任项目重分析/ANA&anti-dsDNA")
####例1
# annotation <- data.frame(group=c(rep('gene',175),rep('pathway',14)))
# rownames(annotation) <- c(rownames(ctrl),rownames(metgene))
group<-c("spearman",'pearson')
for(k in group){
  # data1<-read.xlsx("data.xlsx",sheet = 3,colNames = T,rowNames = T,check.names = F,sep.names = " ")##列是物质，行是样本
  # data1<-as.data.frame(t(data1))  ##data1是列
  # data2<-read.xlsx("data.xlsx",sheet = 3,colNames = T,rowNames = T,check.names = F,sep.names = " ")
  # data2<-as.data.frame(t(data2))
  # annotation<-read.xlsx("data.xlsx",sheet = 4,colNames = T,rowNames = T,check.names = F,sep.names = " ")
  # data2<-data2[rownames(data1),]
  wb <-createWorkbook()
  cor<-matrix(nrow = length(data1),ncol = length(data2))
  cor<-data.frame(cor)
  rownames(cor)<-colnames(data1)
  colnames(cor)<-colnames(data2)
  pvalue<-cor
  for(i in 1:length(data1)){
    for(j in 1:length(data2)){
      x<-as.data.frame(as.data.frame(data1[,i])[!is.na(data1[,i]),])
      y<-as.data.frame(as.data.frame(data2[,j])[!is.na(data2[,j]),])
      if(nrow(x)!=nrow(y)){
        rownames(x)<-rownames(data1)[!is.na(data1[,i])]
        rownames(y)<-rownames(data2)[!is.na(data2[,j])]
        if(nrow(x)<nrow(y)){
        y<-y[rownames(x),]
        x<-as.numeric(x[,1])
      }else{
        x<-x[rownames(y),]
        y<-as.numeric(y[,1])
      }
      }else{
        x<-x[,1]
        y<-y[,1]
      }
      result<-cor.test(x,y,method = k)
      cor[i,j]<-result$estimate
      pvalue[i,j]<-result$p.value
    }
  }
  p<-NULL
  a<-1
  qvalue<-matrix(nrow = nrow(pvalue),ncol = ncol(pvalue))
  qvalue<-data.frame(qvalue)
  rownames(qvalue)<-rownames(pvalue)
  colnames(qvalue)<-colnames(pvalue)
  for(kk in 1:length(pvalue)){
    z<-pvalue[,kk]
    p<-c(p,z)}
  q<-p.adjust(p,method = "BH")
  # q<-NULL
  # prank<-rank(p)
  # for(kkk in 1:length(prank)){
  #   q<-c(q,p[kkk]*(length(p)/prank[kkk]))
  # }
  for(ii in 1:ncol(pvalue)){
    for(jj in 1:nrow(pvalue)){
      qvalue[jj,ii]<-q[a]
      a<-a+1
    }
  }
  X1<-NULL
  X2<-NULL
  Corr<-NULL
  Pvalue<-NULL
  Qvalue<-NULL
  for(m in 1:nrow(pvalue)){
    for(n in 1:ncol(pvalue)){
      if(pvalue[m,n]<1){
        X1<-c(X1,rownames(pvalue)[m])
        X2<-c(X2,colnames(pvalue)[n])
        Corr<-c(Corr,cor[m,n])
        Pvalue<-c(Pvalue,pvalue[m,n])
        Qvalue<-c(Qvalue,qvalue[m,n])
        select<-data.frame(X1,X2,Corr,Pvalue,Qvalue)
      }
    }
  }
  # label<-rownames(annotation)
  # group<-annotation[,1]
  # size<-NULL
  # for(i in label){
  #   if(i %in% colnames(pvalue)){
  #     size<-c(size,length(pvalue[,i][pvalue[,i]<1]))
  #   }else{
  #     size<-c(size,length(unlist(pvalue[i,])[unlist(pvalue[i,])<1]))
  #   }
  # }
  # id<-1:length(label)
  # nodes<-data.frame(id=id,label=label,group=group,size=size)
  # 
  # from<-NULL
  # to<-NULL
  # fromn<-NULL
  # ton<-NULL
  # weight<-NULL
  # for(i in colnames(cor)){
  #   for(j in rownames(cor)){
  #     if(pvalue[j,i]<1){
  #       weight<-c(weight,cor[j,i])
  #       from<-c(from,nodes$id[match(i,nodes$label)])
  #       to<-c(to,nodes$id[match(j,nodes$label)])
  #       fromn<-c(fromn,i)
  #       ton<-c(ton,j)
  #     }
  #   }
  # }
  # edges<-data.frame(from=from,to=to,weight=weight)
  # edges$color<-ifelse(edges$weight<0,"green","purple")
  # edges$weight<-abs(edges$weight)
  # edges<-edges[match(unique(edges$weight),edges$weight),]
  # 
  # edgesc<-data.frame(fromn=fromn,ton=ton,weight=weight)
  # edgesc$color<-ifelse(edgesc$weight<0,"green","purple")
  # edgesc<-edgesc[match(unique(edgesc$weight),edgesc$weight),]
  # 
  addWorksheet(wb,sheetName = "cor",gridLines = F)
  addWorksheet(wb,sheetName = "pvalue",gridLines = F)
  addWorksheet(wb,sheetName = "qvalue",gridLines = F)
  addWorksheet(wb,sheetName = "select",gridLines = F)
  # addWorksheet(wb,sheetName = "nodes",gridLines = F)
  # addWorksheet(wb,sheetName = "edges",gridLines = F)
  # addWorksheet(wb,sheetName = "edgesc",gridLines = F)
  writeDataTable(wb,sheet = 1,cor,rowNames = T)
  writeDataTable(wb,sheet = 2,pvalue,rowNames = T)
  writeDataTable(wb,sheet = 3,qvalue,rowNames = T)
  writeDataTable(wb,sheet = 4,select,rowNames = F)
  # writeDataTable(wb,sheet = 5,nodes,rowNames = F)
  # writeDataTable(wb,sheet = 6,edges,rowNames = F)
  # writeDataTable(wb,sheet = 7,edgesc,rowNames = F)
  saveWorkbook(wb,overwrite = T,paste0(k,".xlsx"))
  limits<-floor(max(max(cor,na.rm = T),1*min(cor,na.rm = T))*100000)/100000
  bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
  c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
  pheatmap(cor,scale = "none",color = c,show_colnames=T,cluster_cols=T,cluster_rows = T,display_numbers = ifelse(qvalue<0.05,"*",""),#annotation_col = annotation,annotation_row = annotation,
           breaks=1*bk,cellwidth = 15, cellheight = 15, fontsize_row  = 12, fontsize_col = 12, border_color = "white",show_rownames = T,fontsize_number = 18,
           filename = paste0(k," 相关性heatmap.png"),annotation_colors = annotation_color)
  pheatmap(cor,scale = "none",color = c,show_colnames=T,cluster_cols=T,cluster_rows = T,display_numbers = ifelse(qvalue<0.05,"*",""),#annotation_col = annotation,annotation_row = annotation,
           breaks=1*bk,cellwidth = 15, cellheight = 15, fontsize_row  = 12, fontsize_col = 12, border_color = "white",show_rownames = T,fontsize_number = 18,
           filename = paste0(k," 相关性heatmap.pdf"),annotation_colors = annotation_color)
}

###选择相关性热图P值
# row<-NULL
# col<-NULL
# for(ii in rownames(cor)){
#   if(any(pvalue[ii,]<0.05)){
#     row<-c(row,ii)
#   }
# }
# for(ii in colnames(cor)){
#   if(any(pvalue[,ii]<0.05)){
#     col<-c(col,ii)
#   }
# }
# cor<-cor[row,col]
# pvalue<-pvalue[row,col]

Enzyme<-jet(5)
names(Enzyme)<-unique(annotation[,1])
annotation_color<-list(Enzyme=Enzyme)
limits<-floor(max(max(cor,na.rm = T),1*min(cor,na.rm = T))*100000)/100000
bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
pheatmap(cor,scale = "none",color = c,show_colnames=T,cluster_cols=T,cluster_rows = T,display_numbers = ifelse(qvalue<0.05,"*",""),#annotation_col = annotation,annotation_row = annotation,
         breaks=1*bk,cellwidth = 15, cellheight = 15, fontsize_row  = 12, fontsize_col = 12, border_color = "white",show_rownames = T,fontsize_number = 18,
         filename = paste0(k," 相关性heatmap.png"),annotation_colors = annotation_color)
pheatmap(cor,scale = "none",color = c,show_colnames=T,cluster_cols=T,cluster_rows = T,display_numbers = ifelse(qvalue<0.05,"*",""),#annotation_col = annotation,annotation_row = annotation,
         breaks=1*bk,cellwidth = 15, cellheight = 15, fontsize_row  = 12, fontsize_col = 12, border_color = "white",show_rownames = T,fontsize_number = 18,
         filename = paste0(k," 相关性heatmap.pdf"),annotation_colors = annotation_color)

####例2
col1 <- colorRampPalette(c("blue","white","red"))
png(paste0(k," 三角相关性heatmap.png"),res = 600,width = 2400,height = 2400)
p<-corrplot::corrplot(as.matrix(cor),method = "square",type = "lower",col = col1(100),tl.pos = "n",order = "hclust",hclust.method = "complete",addrect = 4,rect.col = "red",addgrid.col = "white")
dev.off()
png(paste0(k," 三角相关性heatmap.pdf"),res = 600,width = 2400,height = 2400)
p<-corrplot::corrplot(as.matrix(cor),method = "square",type = "lower",col = col1(100),tl.pos = "n",order = "hclust",hclust.method = "complete",addrect = 4,rect.col = "red",addgrid.col = "white")
dev.off()


###网络图
colors_node<-jet(length(unique(nodes$group)))
names(colors_node)<-unique(nodes$group)
colors_edge<-c("purple","green")
names(colors_edge)<-c("purple","green")
net.tidy <- tbl_graph(
  nodes = nodes, edges = edges, directed = TRUE
)

pdf("stress-q0.05.pdf",width = 30,height = 25)
layout <- create_layout(net.tidy, layout = 'stress')
ggraph(layout) +
# ggraph(net.tidy,layout = "stress")+
  geom_edge_link(aes(edge_width=weight,color=color),show.legend=FALSE) + 
  geom_node_point(aes(fill=group,size=size),shape=21,
                  show.legend = T) +
  scale_edge_width_continuous(range = c(0.5,2)) +
  scale_size(range = c(4, 20))+
  geom_node_text(aes(label = label),repel = TRUE,size=15) +
  scale_fill_manual(values = colors_node)+
  scale_edge_colour_manual(values = colors_edge)+
  theme_graph()
dev.off()

pdf("kk-q0.05.pdf",width = 30,height = 25)
# layout <- create_layout(net.tidy, layout = 'stress')
# ggraph(layout) +
  ggraph(net.tidy,layout = "kk")+
  geom_edge_link(aes(edge_width=weight,color=color),show.legend=FALSE) + 
  geom_node_point(aes(fill=group,size=size),shape=21,
                  show.legend = T) +
  scale_edge_width_continuous(range = c(0.5,2)) +
  scale_size(range = c(4, 20))+
  geom_node_text(aes(label = label),repel = TRUE,size=15) +
  scale_fill_manual(values = colors_node)+
  scale_edge_colour_manual(values = colors_edge)+
  theme_graph()
dev.off()

pdf("centrality-q0.05.pdf",width = 30,height = 25)
# layout <- create_layout(net.tidy, layout = 'stress')
# ggraph(layout) +
ggraph(net.tidy,layout = "centrality",cent=size)+
  geom_edge_link(aes(edge_width=weight,color=color),show.legend=FALSE) + 
  geom_node_point(aes(fill=group,size=size),shape=21,
                  show.legend = T) +
  scale_edge_width_continuous(range = c(0.5,2)) +
  scale_size(range = c(4, 20))+
  geom_node_text(aes(label = label),repel = TRUE,size=15) +
  scale_fill_manual(values = colors_node)+
  scale_edge_colour_manual(values = colors_edge)+
  theme_graph()
dev.off()
  
