library(lattice)
library(latticeExtra)
library(openxlsx)
library(ggplot2)
library(readxl)

####lattice####
# a<-c(1,2,3,4,5)
# b<-c(5,4,3,2,1)
# data<-data.frame(a=a,b=b)
# temp<-data
# hc=hclust(dist(temp)) ###按行聚类
# 
# dd.row=as.dendrogram(hc)###保存行聚类树形
# 
# row.ord=order.dendrogram(dd.row) ###保存行聚类顺序
# 
# hc=hclust(dist(t(temp))) ###按列聚类
# 
# dd.col=as.dendrogram(hc) ###保存列聚类树形
# 
# col.rod=order.dendrogram(dd.col) ###保存列聚类顺序
# 
# temp1=temp[row.ord,] ###只对行聚类（是否对行、列聚类）
# 
# levelplot(t(temp1),aspect="fill",colorkey=list(space="bottom",width=1),xlab="",
#           ylab="",legend=list(right=list(fun=dendrogramGrob,args=list(x=dd.row,rod=row.ord,side='right',size=1)),
#             scales=list(x=list(rot=90))
#             ))###x轴标签旋转60度

orderpro<-read.xlsx("样本顺序.xlsx",sheet = 1,colNames = T,rowNames = T)
ordercli<-read.xlsx("样本顺序.xlsx",sheet = 2,colNames = T,rowNames = T)
data1<-data.frame(sample=ordercli[,"RA"])
for(i in excel_sheets("RA_clinical data_validatoin(1).xlsx")){
  raw<-read.xlsx("RA_clinical data_validatoin(1).xlsx",sheet = i,rowNames = F)
  data<-as.data.frame(raw[,1])
  colnames(data)<-i
  data1<-cbind(data1,data)
}
rownames(data1)<-data1$sample
data1<-data1[,-1]
data1<-data1[orderpro[,"RA"],]
write.xlsx(data1,"RA-clinical.xlsx",rowNames=T,overwrite = T)


####相关性分析####

data2<-read.xlsx("SLE_S-clinical.xlsx",sheet = 1,colNames = T,rowNames = T,check.names = F,sep.names = " ")
# for(ii in 1:ncol(data2)){
#   result<-sapply(as.numeric(data2[,ii]),function(i){(i-mean(as.numeric(data2[,ii])))/sd(as.numeric(data2[,ii]))})
#   data2[,ii]<-result
# }
data1<-read.xlsx("PROTEIN.xlsx",sheet = 1,colNames = T,rowNames = T,check.names = F,sep.names = " ")##列是物质，行是样本
data1<-data1[rownames(data2),]
# for(ii in 1:ncol(data1)){
#   result<-sapply(as.numeric(data1[,ii]),function(i){(i-mean(as.numeric(data1[,ii])))/sd(as.numeric(data1[,ii]))})
#   data1[,ii]<-result
# }
group<-c("spearman",'pearson')
for(k in group){
  # annotation<-read.xlsx("data.xlsx",sheet = 4,colNames = T,rowNames = T,check.names = F,sep.names = " ")

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
  data<-select
  data$cor<-"POS"
  data$cor[data$Corr<0]<-"NEG"
  
  p<-ggplot(data,aes(X1,X2))+
    geom_tile(color="white",aes(fill=-log10(Pvalue)),size=1)+
    scale_fill_gradient(low = "white",high = "black")+
    geom_point(aes(size=abs(Corr),color=cor),shape=16)+
    scale_color_manual(
      name = NULL, values = c(
        'POS' = "#FF7F24", 'NEG' = "#1E90FF"
      ), guide = guide_legend(order = 1, override.aes = list(size = 3))
    )+
    scale_size_continuous(
      name = 'Corr', guide = guide_legend(order = 2), range = c(1.5, 6),
      breaks = c(min(abs(data$Corr)), max(abs(data$Corr))),
      labels = c('0.0', round(max(data$Corr), 1))
    )+
    labs(x = "CellType",y = "Model Gene",color=NULL) + 
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(text=element_text(family="sans"),
          axis.ticks.x = element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(angle = 90,hjust = 1,size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          panel.border = element_rect(fill=NA,color="white",
                                      size=2, linetype="solid"))
  ggsave(
    paste0(k,"_",'all.pdf'), p,
    width = 10, height = 7.5, units = 'in', dpi = 600
  )
  ggsave(
    paste0(k,"_",'all.pdf'), p,
    width = 10, height = 7.5, units = 'in'
  )
  
}





data<-select
data$cor<-"POS"
data$cor[data$Corr<0]<-"NEG"

p<-ggplot(data,aes(X1,X2))+
  geom_tile(color="white",aes(fill=-log10(Pvalue)),size=1)+
  scale_fill_gradient(low = "white",high = "black")+
  geom_point(aes(size=abs(Corr),color=cor),shape=16)+
  scale_color_manual(
    name = NULL, values = c(
      'POS' = "#FF7F24", 'NEG' = "#1E90FF"
    ), guide = guide_legend(order = 1, override.aes = list(size = 3))
  )+
  scale_size_continuous(
    name = 'Corr', guide = guide_legend(order = 2), range = c(1.5, 6),
    breaks = c(min(abs(data$Corr)), max(abs(data$Corr))),
    labels = c('0.0', round(max(data$Corr), 1))
  )+
  labs(x = "Protein",y = "Clinical",color=NULL) + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(family="sans"),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill=NA,color="white",
                                    size=2, linetype="solid"))
ggsave(
  'SLE_S.jpg', p,
  width = 10, height = 7.5, units = 'in', dpi = 600
)
ggsave(
  'SLE_S.pdf', p,
  width = 10, height = 7.5, units = 'in'
)





