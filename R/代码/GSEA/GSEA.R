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

group1<-list.files()[6:8]
# group1<-c("SLE VS HC_RFA.xlsx","RA VS HC_RFA.xlsx")
for(i in group1){
  data<-read.xlsx(i,sheet = 1,rowNames = T,sep.names = " ",check.names = F)
  temp1<-unlist(strsplit(i,"_RFA.xlsx"))
  if (!dir.exists(temp1)){dir.create(temp1)}
  data1<-data[,grep(strsplit(temp1," VS ")[[1]][1],colnames(data))]
  data2<-data[,grep(strsplit(temp1," VS ")[[1]][2],colnames(data))]
  data1$mean<-sapply(rownames(data1),function(x){mean(unlist(data1[x,]))})
  data2$mean<-sapply(rownames(data2),function(x){mean(unlist(data2[x,]))})
  gene<-data.frame(SYMBOL=rownames(data),log2fc=log2(data1$mean/data2$mean))
  genetran<-bitr(gene[,1],fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  genetran<-dplyr::distinct(genetran,SYMBOL,.keep_all=TRUE)
  genedt<-merge(gene,genetran,by = "SYMBOL")
  write.xlsx(genedt,paste0(temp1,"/","symbol-entrezid-logfc.xlsx"))
  genelist<-genedt$log2fc
  names(genelist)<-genedt$ENTREZID
  genelist<-sort(genelist,decreasing = T)
  # group2<-list.files()[1:4]
  # group2<-c("c5.go.cc.v7.4.entrez.txt","c5.go.bp.v7.4.entrez.txt","c5.go.mf.v7.4.entrez.txt","c2.cp.kegg.v7.4.entrez.txt")
  group2<-c("c2.cp.kegg.v7.4.entrez.txt")
  for(j in group2){
    temp2<-strsplit(j,"\\.")[[1]][3]
    if (!dir.exists(paste0(temp1,"/",temp2))){dir.create(paste0(temp1,"/",temp2))}
    gmt<-read.gmt(j)
    gsea<-GSEA(genelist,TERM2GENE = gmt,minGSSize = 1,maxGSSize = 500,pvalueCutoff = 1)
    gsearesult<-gsea@result
    gsearesult$count<-sapply(gsearesult$core_enrichment,function(x){length(unlist(str_split(x,"/")))})
    gsearesult$color<-ifelse(gsearesult$NES<0,"suppressed(p<0.05)","activated(p<0.05)")
    gsearesult$color[gsearesult$color=="suppressed(p<0.05)"&gsearesult$p.adjust<0.05]<-"suppressed(p.adjust<0.05)"
    gsearesult$color[gsearesult$color=="activated(p<0.05)"&gsearesult$p.adjust<0.05]<-"activated(p.adjust<0.05)"
    gsearesult$color[gsearesult$pvalue>0.05]<-"not significant"
    write.xlsx(gsearesult,paste0(temp1,"/",temp2,"/",j,"-result.xlsx"))
    save(gsea,file = paste0(temp1,"/",temp2,"/",j,"-result.Rdata"))
    temp3<-"显著性前50富集曲线"
    if (!dir.exists(paste0(temp1,"/",temp2,"/",temp3))){dir.create(paste0(temp1,"/",temp2,"/",temp3))}
    if(length(gsea@result$Description)>50){
      kvalue<-1:50
    }else{
      kvalue<-1:length(gsea@result$Description)
    }
    for(k in kvalue){
      
      png(paste0(temp1,"/",temp2,"/",temp3,"/",gsea@result$Description[k],".png"),res = 600,width = 7500,height = 6000)
      p<-gseaplot2(gsea,k,color="red",pvalue_table = T,base_size = 24)
      print(p)
      dev.off()
      
      pdf(paste0(temp1,"/",temp2,"/",temp3,"/",gsea@result$Description[k],".pdf"),width = 15,height = 12)
      p<-gseaplot2(gsea,k,color="red",pvalue_table = T,base_size = 24)
      print(p)
      dev.off()
      
    }
    if(nrow(gsearesult[gsearesult$pvalue<1,])>20){
      pnumber<-1:20
    }else{
      pnumber<-1:nrow(gsearesult[gsearesult$pvalue<1,])
    }
    pathwayplot<-gsearesult[pnumber,]
    pathwayplot<-pathwayplot[order(pathwayplot$NES),]
    pathwayplot$ID<-factor(pathwayplot$ID,levels = pathwayplot$ID)
    color<-c("#FF3030","#FF303050","#4876FF","#4876FF50","#696969")
    names(color)<-c("activated(p.adjust<0.05)","activated(p<0.05)","suppressed(p.adjust<0.05)","suppressed(p<0.05)","not significant")

    # color<-c("#FF303050","#696969")
    # names(color)<-c("activated(p<0.05)","not significant")
    
    p<-ggplot(pathwayplot,aes(NES,ID))+
      geom_point(
        aes(NES, ID, size = count,col = color)
      )+
      scale_color_manual(values = color)+
      scale_size_continuous(
        name = 'Hits Number', range = c(4, 12),
        breaks = c(min(pathwayplot$count), max(pathwayplot$count)),
        labels = c(min(pathwayplot$count), round(max(pathwayplot$count), 4)),
        guide = guide_legend(title.hjust = 0.5)
      )+
      theme_bw()+
      theme(
        panel.grid.major = element_line(
          color = 'grey', linetype = 'dashed'
        ),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = c(.95, .05),
        legend.justification = c(.95, .05)
      )+
      labs(
        x = "NES",
        y = "Pathway",
        size = "Hits Number"
      )
    ggsave(paste0(temp1,"/",temp2,"/","pathwayplot.png"),p,dpi = 600,width = 0.75*nrow(pathwayplot),height = 9,limitsize = F,units = "in")
    ggsave(paste0(temp1,"/",temp2,"/","pathwayplot.pdf"),p,width = 0.75*nrow(pathwayplot),height = 9,limitsize = F,units = "in")
    pdf(paste0(temp1,"/",temp2,"/","波浪图.pdf"),width=10,height=12)
    p<-ridgeplot(gsea,10)
    print(p)
    dev.off()
    pdf(paste0(temp1,"/",temp2,"/","气泡图-激活抑制.pdf"),width=10,height=9)
    p<-dotplot(gsea,split=".sign")+facet_grid(~.sign) 
    print(p)
    dev.off()
    pdf(paste0(temp1,"/",temp2,"/","气泡图.pdf"),width=10,height=12)
    p<-dotplot(gsea)
    print(p)
    dev.off()
    pdf(paste0(temp1,"/",temp2,"/","富集曲线前10NES.pdf"),width=15,height=12)
    p<-gseaplot2(gsea,1:10,color="red",pvalue_table = T,base_size = 24)
    print(p)
    dev.off()
  }
}




