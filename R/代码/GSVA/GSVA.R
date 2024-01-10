library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(GSVAdata)
library(DOSE)
library(stringr)
library(openxlsx)
library(ggplot2)
library(ggridges)
library(dplyr)
library(pals)
library(pheatmap)
library(tcltk)
library(affy)
library(MKomics)

setwd("D:/马驰宇/181/戴主任/代码/GSVA/entrezid")


###用limma包差异分析
# group1<-list.files()[5:7]
group1<-c("DE protein.xlsx")
for(i in group1){
  data<-read.xlsx(i,sheet = 3,rowNames = T,sep.names = " ",check.names = F)
  annotation<-read.xlsx(i,sheet = 7,rowNames = T,sep.names = " ",check.names = F)
  temp1<-unlist(strsplit(i,".xlsx"))
  if (!dir.exists(temp1)){dir.create(temp1)}
  data$SYMBOL<-rownames(data)
  genetran<-bitr(rownames(data),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  genetran<-dplyr::distinct(genetran,SYMBOL,.keep_all=TRUE)
  # expr<-merge(data,genetran,by = "SYMBOL")
  expr<-inner_join(data,genetran,by=c("SYMBOL"="SYMBOL"))
  rownames(expr)<-expr$ENTREZID
  dataexpr<-expr[,1:(ncol(expr)-2)]

  group2<-list.files()[1:4]
  # group2<-c("c5.go.cc.v7.4.entrez.txt","c5.go.bp.v7.4.entrez.txt","c5.go.mf.v7.4.entrez.txt")
  for(j in group2){
    temp2<-strsplit(j,"\\.")[[1]][3]
    if (!dir.exists(paste0(temp1,"/",temp2))){dir.create(paste0(temp1,"/",temp2))}
    gmt<-read.gmt(j)
    gmt_list<-split(gmt$gene,gmt$term)
    gsva<-gsva(as.matrix(dataexpr),gmt_list,kcdf="Gaussian",method = "gsva",min.sz =1)
    gsva<-as.data.frame(gsva)

    design <- model.matrix(~0 + factor(c(rep("Control",15),rep("RA",16),rep("SLE",21))))
    colnames(design) <- c("Control", "RA","SLE")
    row.names(design)<-colnames(gsva)

    temp3<-"联合gsva联合差异"
    padjust<-NULL
    gsvaresult<-gsva
    if (!dir.exists(paste0(temp1,"/",temp2,"/",temp3))){dir.create(paste0(temp1,"/",temp2,"/",temp3))}
    for(k in c("SLE-RA","SLE-Control","RA-Control")){
      contrast.matrix<-makeContrasts(k,levels=design)
      fit1 <- lmFit(gsva, design)
      fit2 <- contrasts.fit(fit1, contrast.matrix)
      fit3 <- eBayes(fit2)
      #这边是总的
      allGeneSets <- topTable(fit3, coef=1, number=Inf,adjust.method = "BH",sort.by = NULL)

      allGeneSets$type<-"not significant"
      allGeneSets$type[allGeneSets$adj.P.Val<0.001&allGeneSets$logFC<0]<-"down"
      allGeneSets$type[allGeneSets$adj.P.Val<0.001&allGeneSets$logFC>0]<-"up"
      allGeneSets$type<-factor(allGeneSets$type,levels = c('down', 'not significant', 'up'))
      p<-ggplot()+
        geom_blank(data = allGeneSets,aes(-logFC, -log10(adj.P.Val)))+
        geom_hline(
          yintercept = c(-log10(0.001)),
          color = "grey50",linetype = "dashed"
        )+
        geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
        geom_point(
          data = allGeneSets[allGeneSets$type=='not significant',],
          aes(logFC,-log10(adj.P.Val), color = type),size = 2.4
        )+
        geom_point(
          data = allGeneSets[allGeneSets$type!='not significant',],
          aes(logFC,-log10(adj.P.Val), color = type),size = 2.4
        )+
        scale_color_manual(
          name = "Status", values = c(
            'down' = '#619cffa0', 'not significant' = '#b3b3b350',
            'up' = '#f8766da0'
          ), guide = guide_legend(order = 1, override.aes = list(size = 3))
        )+
        # scale_size_continuous(
        #   name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6),
        #   breaks = c(min(VolcanoPlot$vip), max(VolcanoPlot$vip)),
        #   labels = c('0.0', round(max(VolcanoPlot$vip), 1))
        # )+
        theme_bw()+
        theme(
          # legend.key = element_blank()
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black",size = 1),
          legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
          legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
          legend.key.size = unit(0.8,"cm"),
          axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
          axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
        )+
        labs(
          x = expression(paste(plain("log"["2"]), ' Fold Change')),
          y = expression(paste(plain("-log"["10"]), ' ', italic('P.adjust'), '-value'))
        )
      ggsave(
        paste0(temp1,"/",temp2,"/",temp3,"/",k,'-volcano plot.jpg'), p,
        width = 9, height = 7.5, units = 'in', dpi = 600
      )
      ggsave(
        paste0(temp1,"/",temp2,"/",temp3,"/",k,'-volcano plot.pdf'), p,
        width = 9, height = 7.5, units = 'in', dpi = 600
      )
      padjust<-c(padjust,rownames(allGeneSets)[allGeneSets$adj.P.Val<0.001&abs(allGeneSets$logFC)>log2(1)])
      colnames(allGeneSets)<-paste0(k,"-",colnames(allGeneSets))
      allGeneSets<-allGeneSets[rownames(gsva),]
      gsvaresult<-cbind(gsvaresult,allGeneSets)
    }
    padjust<-unique(padjust)
    write.xlsx(gsvaresult,paste0(temp1,"/",temp2,"/",temp3,"/","gsvaresult.xlsx"),rowNames=T)
    gsvamap<-gsva[,colnames(dataexpr)]
    gsvamapp<-gsva[padjust,colnames(dataexpr)]
    Sample<-c("#EA70AD","#6BB5DC","#A25BC3")
    names(Sample)<-unique(annotation[,1])
    annotationcolor<-list(Sample=Sample)
    ##热图每个格子的颜色
    limits<-floor(max(max(gsvamap,na.rm = T),1*min(gsvamap,na.rm = T))*100000)/100000
    bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
    c<-c(colorRampPalette(colors = c("#11675F","white"))(length(bk)/2),colorRampPalette(colors = c("white","#8E4E17"))(length(bk)/2))
    ##pheatmap函数绘制热图
    pheatmap(gsvamap,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamap), cellheight = 400/nrow(gsvamap), fontsize_row  = 200/ncol(gsvamap), fontsize_col = 200/ncol(gsvamap), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-cluster-gsvamap.png"))
    pheatmap(gsvamap,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamap), cellheight = 400/nrow(gsvamap), fontsize_row  = 200/ncol(gsvamap), fontsize_col = 200/ncol(gsvamap), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-cluster-gsvamap.pdf"))
    pheatmap(gsvamap,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamap), cellheight = 400/nrow(gsvamap), fontsize_row  = 200/ncol(gsvamap), fontsize_col = 200/ncol(gsvamap), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-nocluster-gsvamap.png"))
    pheatmap(gsvamap,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamap), cellheight = 400/nrow(gsvamap), fontsize_row  = 200/ncol(gsvamap), fontsize_col = 200/ncol(gsvamap), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-nocluster-gsvamap.pdf"))
    ##热图每个格子的颜色
    limits<-floor(max(max(gsvamapp,na.rm = T),1*min(gsvamapp,na.rm = T))*100000)/100000
    bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
    c<-c(colorRampPalette(colors = c("#11675F","white"))(length(bk)/2),colorRampPalette(colors = c("white","#8E4E17"))(length(bk)/2))
    ##pheatmap函数绘制热图
    pheatmap(gsvamapp,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamapp), cellheight = 400/nrow(gsvamapp), fontsize_row  = 200/ncol(gsvamapp), fontsize_col = 200/ncol(gsvamapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-cluster-gsvamapp.png"))
    pheatmap(gsvamapp,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamapp), cellheight = 400/nrow(gsvamapp), fontsize_row  = 200/ncol(gsvamapp), fontsize_col = 200/ncol(gsvamapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-cluster-gsvamapp.pdf"))
    pheatmap(gsvamapp,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamapp), cellheight = 400/nrow(gsvamapp), fontsize_row  = 200/ncol(gsvamapp), fontsize_col = 200/ncol(gsvamapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-nocluster-gsvamapp.png"))
    pheatmap(gsvamapp,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsvamapp), cellheight = 400/nrow(gsvamapp), fontsize_row  = 200/ncol(gsvamapp), fontsize_col = 200/ncol(gsvamapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp3,"/","-nocluster-gsvamapp.pdf"))

    temp4<-"联合gsva独立差异"
    if (!dir.exists(paste0(temp1,"/",temp2,"/",temp4))){dir.create(paste0(temp1,"/",temp2,"/",temp4))}
    padjust<-NULL
    gsva2result<-gsva
    for(k in c("SLE-RA","SLE-Control","RA-Control")){
      gsva2<-gsva[,grep(str_replace(k,"-","|"),colnames(gsva))]
      design <- model.matrix(~0 + factor(annotation[colnames(gsva2),1]))
      colnames(design) <- unique(annotation[colnames(gsva2),1])
      row.names(design)<-colnames(gsva2)
      contrast.matrix<-makeContrasts(k,levels=design)
      fit1 <- lmFit(gsva2, design)
      fit2 <- contrasts.fit(fit1, contrast.matrix)
      fit3 <- eBayes(fit2)
      #这边是总的
      allGeneSets <- topTable(fit3, coef=1, number=Inf,adjust.method = "BH",sort.by = NULL)

      allGeneSets$type<-"not significant"
      allGeneSets$type[allGeneSets$adj.P.Val<0.001&allGeneSets$logFC<0]<-"down"
      allGeneSets$type[allGeneSets$adj.P.Val<0.001&allGeneSets$logFC>0]<-"up"
      allGeneSets$type<-factor(allGeneSets$type,levels = c('down', 'not significant', 'up'))
      p<-ggplot()+
        geom_blank(data = allGeneSets,aes(-logFC, -log10(adj.P.Val)))+
        geom_hline(
          yintercept = c(-log10(0.001)),
          color = "grey50",linetype = "dashed"
        )+
        geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
        geom_point(
          data = allGeneSets[allGeneSets$type=='not significant',],
          aes(logFC,-log10(adj.P.Val), color = type),size = 2.4
        )+
        geom_point(
          data = allGeneSets[allGeneSets$type!='not significant',],
          aes(logFC,-log10(adj.P.Val), color = type),size = 2.4
        )+
        scale_color_manual(
          name = "Status", values = c(
            'down' = '#619cffa0', 'not significant' = '#b3b3b350',
            'up' = '#f8766da0'
          ), guide = guide_legend(order = 1, override.aes = list(size = 3))
        )+
        # scale_size_continuous(
        #   name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6),
        #   breaks = c(min(VolcanoPlot$vip), max(VolcanoPlot$vip)),
        #   labels = c('0.0', round(max(VolcanoPlot$vip), 1))
        # )+
        theme_bw()+
        theme(
          # legend.key = element_blank()
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black",size = 1),
          legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
          legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
          legend.key.size = unit(0.8,"cm"),
          axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
          axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
        )+
        labs(
          x = expression(paste(plain("log"["2"]), ' Fold Change')),
          y = expression(paste(plain("-log"["10"]), ' ', italic('P.adjust'), '-value'))
        )
      ggsave(
        paste0(temp1,"/",temp2,"/",temp4,"/",k,'-volcano plot.jpg'), p,
        width = 9, height = 7.5, units = 'in', dpi = 600
      )
      ggsave(
        paste0(temp1,"/",temp2,"/",temp4,"/",k,'-volcano plot.pdf'), p,
        width = 9, height = 7.5, units = 'in', dpi = 600
      )
      padjust<-c(padjust,rownames(allGeneSets)[allGeneSets$adj.P.Val<0.001&abs(allGeneSets$logFC)>log2(1)])
      colnames(allGeneSets)<-paste0(k,"-",colnames(allGeneSets))
      allGeneSets<-allGeneSets[rownames(gsva2),]
      gsva2result<-cbind(gsva2result,allGeneSets)
    }
    padjust<-unique(padjust)
    write.xlsx(gsva2result,paste0(temp1,"/",temp2,"/",temp4,"/","gsvaresult.xlsx"),rowNames=T)
    gsva2map<-gsva2result[,colnames(dataexpr)]
    gsva2mapp<-gsva2result[padjust,colnames(dataexpr)]
    Sample<-c("#EA70AD","#6BB5DC","#A25BC3")
    names(Sample)<-unique(annotation[,1])
    annotationcolor<-list(Sample=Sample)
    ##热图每个格子的颜色
    limits<-floor(max(max(gsva2map,na.rm = T),1*min(gsva2map,na.rm = T))*100000)/100000
    bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
    c<-c(colorRampPalette(colors = c("#11675F","white"))(length(bk)/2),colorRampPalette(colors = c("white","#8E4E17"))(length(bk)/2))
    ##pheatmap函数绘制热图
    pheatmap(gsva2map,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2map), cellheight = 400/nrow(gsva2map), fontsize_row  = 200/ncol(gsva2map), fontsize_col = 200/ncol(gsva2map), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-cluster-gsva2map.png"))
    pheatmap(gsva2map,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2map), cellheight = 400/nrow(gsva2map), fontsize_row  = 200/ncol(gsva2map), fontsize_col = 200/ncol(gsva2map), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-cluster-gsva2map.pdf"))
    pheatmap(gsva2map,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2map), cellheight = 400/nrow(gsva2map), fontsize_row  = 200/ncol(gsva2map), fontsize_col = 200/ncol(gsva2map), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-nocluster-gsva2map.png"))
    pheatmap(gsva2map,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2map), cellheight = 400/nrow(gsva2map), fontsize_row  = 200/ncol(gsva2map), fontsize_col = 200/ncol(gsva2map), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-nocluster-gsva2map.pdf"))
    ##热图每个格子的颜色
    limits<-floor(max(max(gsva2mapp,na.rm = T),1*min(gsva2mapp,na.rm = T))*100000)/100000
    bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
    c<-c(colorRampPalette(colors = c("#11675F","white"))(length(bk)/2),colorRampPalette(colors = c("white","#8E4E17"))(length(bk)/2))
    ##pheatmap函数绘制热图
    pheatmap(gsva2mapp,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2mapp), cellheight = 400/nrow(gsva2mapp), fontsize_row  = 200/ncol(gsva2mapp), fontsize_col = 200/ncol(gsva2mapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-cluster-gsva2mapp.png"))
    pheatmap(gsva2mapp,scale = "none",color = c,cluster_cols=T,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2mapp), cellheight = 400/nrow(gsva2mapp), fontsize_row  = 200/ncol(gsva2mapp), fontsize_col = 200/ncol(gsva2mapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-cluster-gsva2mapp.pdf"))
    pheatmap(gsva2mapp,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2mapp), cellheight = 400/nrow(gsva2mapp), fontsize_row  = 200/ncol(gsva2mapp), fontsize_col = 200/ncol(gsva2mapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-nocluster-gsva2mapp.png"))
    pheatmap(gsva2mapp,scale = "none",color = c,cluster_cols=F,cluster_rows = T,annotation_col = annotation,show_rownames = F,
             breaks=1*bk,cellwidth = 200/ncol(gsva2mapp), cellheight = 400/nrow(gsva2mapp), fontsize_row  = 200/ncol(gsva2mapp), fontsize_col = 200/ncol(gsva2mapp), border=FALSE,annotation_colors = annotationcolor,
             filename = paste0(temp1,"/",temp2,"/",temp4,"/","-nocluster-gsva2mapp.pdf"))

    temp5<-"独立gsva独立差异"
    if (!dir.exists(paste0(temp1,"/",temp2,"/",temp5))){dir.create(paste0(temp1,"/",temp2,"/",temp5))}
    padjust<-NULL
    for(k in c("SLE-RA","SLE-Control","RA-Control")){
      dataexpr3<-dataexpr[,grep(str_replace(k,"-","|"),colnames(dataexpr))]
      gsva3<-gsva(as.matrix(dataexpr3),gmt_list,kcdf="Gaussian",method = "gsva")
      gsva3<-as.data.frame(gsva3)
      design <- model.matrix(~0 + factor(annotation[colnames(gsva3),1]))
      colnames(design) <- unique(annotation[colnames(gsva3),1])
      row.names(design)<-colnames(gsva3)
      contrast.matrix<-makeContrasts(k,levels=design)
      fit1 <- lmFit(gsva3, design)
      fit2 <- contrasts.fit(fit1, contrast.matrix)
      fit3 <- eBayes(fit2)
      #这边是总的
      allGeneSets <- topTable(fit3, coef=1, number=Inf,adjust.method = "BH",sort.by = NULL)

      allGeneSets$type<-"not significant"
      allGeneSets$type[allGeneSets$adj.P.Val<0.05&allGeneSets$logFC<0]<-"down"
      allGeneSets$type[allGeneSets$adj.P.Val<0.05&allGeneSets$logFC>0]<-"up"
      allGeneSets$type<-factor(allGeneSets$type,levels = c('down', 'not significant', 'up'))
      p<-ggplot()+
        geom_blank(data = allGeneSets,aes(-logFC, -log10(adj.P.Val)))+
        geom_hline(
          yintercept = c(-log10(0.05)),
          color = "grey50",linetype = "dashed"
        )+
        geom_vline(xintercept = 0, color = 'grey50', linetype = 'dashed') +
        geom_point(
          data = allGeneSets[allGeneSets$type=='not significant',],
          aes(logFC,-log10(adj.P.Val), color = type),size = 2.4
        )+
        geom_point(
          data = allGeneSets[allGeneSets$type!='not significant',],
          aes(logFC,-log10(adj.P.Val), color = type),size = 2.4
        )+
        scale_color_manual(
          name = "Status", values = c(
            'down' = '#619cffa0', 'not significant' = '#b3b3b350',
            'up' = '#f8766da0'
          ), guide = guide_legend(order = 1, override.aes = list(size = 3))
        )+
        # scale_size_continuous(
        #   name = 'VIP', guide = guide_legend(order = 2), range = c(2, 6),
        #   breaks = c(min(VolcanoPlot$vip), max(VolcanoPlot$vip)),
        #   labels = c('0.0', round(max(VolcanoPlot$vip), 1))
        # )+
        theme_bw()+
        theme(
          # legend.key = element_blank()
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black",size = 1),
          legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
          legend.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
          legend.key.size = unit(0.8,"cm"),
          axis.title = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
          axis.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5)
        )+
        labs(
          x = expression(paste(plain("log"["2"]), ' Fold Change')),
          y = expression(paste(plain("-log"["10"]), ' ', italic('P.adjust'), '-value'))
        )
      ggsave(
        paste0(temp1,"/",temp2,"/",temp5,"/",k,'-volcano plot.jpg'), p,
        width = 9, height = 7.5, units = 'in', dpi = 600
      )
      ggsave(
        paste0(temp1,"/",temp2,"/",temp5,"/",k,'-volcano plot.pdf'), p,
        width = 9, height = 7.5, units = 'in', dpi = 600
      )
      padjust<-c(padjust,rownames(allGeneSets)[allGeneSets$adj.P.Val<0.001&abs(allGeneSets$logFC)>log2(1)])
      colnames(allGeneSets)<-paste0(k,"-",colnames(allGeneSets))
      allGeneSets<-allGeneSets[rownames(gsva3),]
      gsva3<-cbind(gsva3,allGeneSets)
      write.xlsx(gsva3,paste0(temp1,"/",temp2,"/",temp5,"/",k,"gsvaresult.xlsx"),rowNames=T)
    }
    padjust<-unique(padjust)
    write.xlsx(gsva3[padjust,],paste0(temp1,"/",temp2,"/",temp5,"/","gsvaresultp.xlsx"),rowNames=T)
  }
}


###双样本T检验和moderated t test
group1<-c("profile10_sample expression.xlsx")
for(i in group1){
  data<-read.xlsx(i,sheet = 1,rowNames = T,sep.names = " ",check.names = F)
  annotation<-read.xlsx(i,sheet = 2,rowNames = T,sep.names = " ",check.names = F)
  temp1<-unlist(strsplit(i,".xlsx"))
  if (!dir.exists(temp1)){dir.create(temp1)}
  data$SYMBOL<-rownames(data)
  genetran<-bitr(rownames(data),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  genetran<-dplyr::distinct(genetran,SYMBOL,.keep_all=TRUE)
  # expr<-merge(data,genetran,by = "SYMBOL")
  expr<-inner_join(data,genetran,by=c("SYMBOL"="SYMBOL"))
  rownames(expr)<-expr$ENTREZID
  dataexpr<-expr[,1:(ncol(expr)-2)]
  
  # group2<-list.files()[1:4]
  # group2<-c("c5.go.cc.v7.4.entrez.txt","c5.go.bp.v7.4.entrez.txt","c5.go.mf.v7.4.entrez.txt")
  group2<-c("c2.cp.kegg.v7.4.entrez.txt")
  for(j in group2){
    temp2<-strsplit(j,"\\.")[[1]][3]
    if (!dir.exists(paste0(temp1,"/",temp2))){dir.create(paste0(temp1,"/",temp2))}
    gmt<-read.gmt(j)
    gmt_list<-split(gmt$gene,gmt$term)
    gsva<-gsva(as.matrix(dataexpr),gmt_list,kcdf="Gaussian",method = "gsva",min.sz=1.1)
    gsva<-as.data.frame(gsva)
    pvalue<-NULL
    for(k in 1:nrow(gsva)){
      var<-var.test(unlist(gsva[k,grep("Control",colnames(gsva))],use.names = F),unlist(gsva[k,grep("SLEA",colnames(gsva))],use.names = F))$p.value > 0.05
      p<-t.test(unlist(gsva[k,grep("Control",colnames(gsva))],use.names = F), unlist(gsva[k,grep("SLEA",colnames(gsva))],use.names = F), var.equal = var)$p.value
      pvalue<-c(pvalue,p)
    }
    qvalue<-p.adjust(pvalue,method = "BH")
    ttest<-data.frame("C-SA-twosample-p"=pvalue,"C-SA-twosample-q"=qvalue)
    gsva<-cbind(gsva,ttest)
    pvalue<-NULL
    for(k in 1:nrow(gsva)){
      var<-var.test(unlist(gsva[k,grep("Control",colnames(gsva))],use.names = F),unlist(gsva[k,grep("SLES",colnames(gsva))],use.names = F))$p.value > 0.05
      p<-t.test(unlist(gsva[k,grep("Control",colnames(gsva))],use.names = F), unlist(gsva[k,grep("SLES",colnames(gsva))],use.names = F), var.equal = var)$p.value
      pvalue<-c(pvalue,p)
    }
    qvalue<-p.adjust(pvalue,method = "BH")
    ttest<-data.frame("C-SS-twosample-p"=pvalue,"C-SS-twosample-q"=qvalue)
    gsva<-cbind(gsva,ttest)
    pvalue<-NULL
    for(k in 1:nrow(gsva)){
      var<-var.test(unlist(gsva[k,grep("SLEA",colnames(gsva))],use.names = F),unlist(gsva[k,grep("SLES",colnames(gsva))],use.names = F))$p.value > 0.05
      p<-t.test(unlist(gsva[k,grep("SLEA",colnames(gsva))],use.names = F), unlist(gsva[k,grep("SLES",colnames(gsva))],use.names = F), var.equal = var)$p.value
      pvalue<-c(pvalue,p)
    }
    qvalue<-p.adjust(pvalue,method = "BH")
    ttest<-data.frame("SA-SS-twosample-p"=pvalue,"SA-SS-twosample-q"=qvalue)
    gsva<-cbind(gsva,ttest)
    ###moderated t test
    padjust<-NULL
    gsva2result<-gsva
    for(k in c("SLEA-SLES","SLES-Control","SLEA-Control")){
      gsva2<-gsva[,grep(str_replace(k,"-","|"),colnames(gsva))]
      design <- model.matrix(~0 + factor(annotation[colnames(gsva2),1]))
      colnames(design) <- unique(annotation[colnames(gsva2),1])
      row.names(design)<-colnames(gsva2)
      contrast.matrix<-makeContrasts(k,levels=design)
      fit1 <- lmFit(gsva2, design)
      fit2 <- contrasts.fit(fit1, contrast.matrix)
      fit3 <- eBayes(fit2)
      #这边是总的
      allGeneSets <- topTable(fit3, coef=1, number=Inf,adjust.method = "BH",sort.by = NULL)
      colnames(allGeneSets)<-paste0(k,"-",colnames(allGeneSets))
      allGeneSets<-allGeneSets[rownames(gsva2),]
      gsva2result<-cbind(gsva2result,allGeneSets)
    }
    write.xlsx(gsva2result,paste0(temp1,"/",temp2,"/","gsvaresult.xlsx"),rowNames=T)
  }
}


##另一个包
# set.seed(123)
# X <- rbind(matrix(rnorm(5*20), nrow = 5, ncol = 20),
#            matrix(rnorm(5*20, mean = 1), nrow = 5, ncol = 20))
# g2 <- factor(c(rep("group 1", 10), rep("group 2", 10)))
# mod.t.test(X, group = g2)
# 
# 
# data<-read.xlsx("gsvaresult.xlsx",rowNames = T)
# dat<-as.matrix(data[,grep("Control|SLE",colnames(data)[1:52])])
# group<-factor(c(rep("Control", 15), rep("SLE", 21)))
# result<-mod.t.test(dat, group = group)
