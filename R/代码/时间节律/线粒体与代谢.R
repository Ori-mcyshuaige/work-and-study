library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
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
library(ggpubr)
library(ReactomePA)
library(reactome.db)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("reactome.db")
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)

GO_cc <- read.csv('09.线粒体节律与代谢/线粒体基因go_cc.csv')
pdata <- list(GO_CC=GO_CC)
for(k in names(pdata)){
  data <- pdata[[k]]
  p <- ggplot(data,aes(x=Group,y=Description,fill=p.adjust,
                       color=p.adjust,size=Count))+
    geom_point(shape=16)+
    # scale_x_discrete(limits = c('cluster1','cluster2','cluster3','cluster4','cluster5','cluster6'))+
    # scale_y_discrete(labels=function(x) str_wrap(x, width=80))+
    scale_color_gradient(
      name = "p.adjust", 
      low = 'red', high = 'blue', 
      guide = guide_colorbar(
        title.position = 'top', title.hjust = 0.5, 
        direction = 'horizontal'
      ))+
    scale_fill_gradient(
      name = "p.adjust", 
      low = 'red', high = 'blue', 
      guide = guide_colorbar(
        title.position = 'top', title.hjust = 0.5, 
        direction = 'horizontal'
      ))+
    scale_size_continuous(
      name = 'Count', guide = guide_legend(order = 2), range = c(3, 9),
      breaks = c(min(abs(data$Count)), max(abs(data$Count))),
      labels = c(round(min(data$Count), 1), round(max(data$Count), 1))
    )+
    theme_bw()+
    theme(#panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      text=element_text(family="sans"),
      axis.ticks.x = element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.x = element_text(angle = 90,hjust = 1,size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      panel.border = element_rect(fill=NA,color="black",
                                  size=1, linetype="solid"))+
    # facet_grid(~data$Group,scales= "free",space= "free")+
    labs(x = NULL,y = k,color=NULL) 
  
  ggsave(
    paste0(k,"_",'all.pdf'), p,
    width = 10, height = 8, units = 'in', dpi = 600
  )
  ggsave(
    paste0(k,"_",'all.png'), p,
    width = 10, height = 8, units = 'in'
  )
}
mig <- sapply(GO_CC[GO_CC$Group=='Dm',]$geneID,function(x){strsplit(x,'/')})
names(mig) <- GO_CC[GO_CC$Group=='Dm',]$Description
mig <- unique(unlist(mig))
hg <- bitr(mig,fromType="ENTREZID",toType=c("ENTREZID","ENSEMBL","SYMBOL"),
           OrgDb="org.Mm.eg.db")## 鼠 org.Mm.eg.db   人org.Hs.eg.db
head(hg)
hg <- hg[match(unique(hg$SYMBOL),hg$SYMBOL),]
write.csv(hg,'gene.csv')
allrhy <- read.csv('01.节律与差异节律/keeplogcpm.csv',row.names = 1)
ctrl <- allrhy[hg$SYMBOL,grep('Ctrl',colnames(allrhy))]
dm <- allrhy[hg$SYMBOL,grep('DM',colnames(allrhy))]
ssgsea <- read.csv('07.ssGSEA代谢评分/ssGSEA.metscore.csv',row.names = 1)
gsva <- read.csv('04.GSVA/共有节律基因KEGG/gsva_res.csv',row.names = 1)
met <- rbind(ssgsea,gsva[c(grep('RHYTHM',rownames(gsva)),grep('FATTY',rownames(gsva))),colnames(ssgsea)])
metgene <- c('Cd73','Nadk','Nadsyn1','Nampt','Naprt','Nmnat3','Nnt','Nrk1',
             'Paps','Pnp','Qprt','Sirts','Acadl','Acsl4','Cact','Cd36','Cpt1a','Cpt2','Crat','Fatp2')
metgene <- allrhy[metgene,]
metgene <- metgene[!is.na(metgene$Ctrl_ZT0_1),]
data1 <- as.data.frame(t(cbind(ctrl,dm)))#cbind(ctrl,dm)
data2 <- as.data.frame(t(metgene))


#######
cor <- read.xlsx('09.线粒体节律与代谢/代谢基因相关性/DM/pearson.xlsx',sheet = 1,rowNames = T,colNames = T,check.names = F)
qvalue <- read.xlsx('09.线粒体节律与代谢/代谢基因相关性/DM/pearson.xlsx',sheet = 3,rowNames = T,colNames = T,check.names = F)
sel <- read.xlsx('09.线粒体节律与代谢/代谢基因相关性/DM/pearson.xlsx',sheet = 4,rowNames = F,colNames = T,check.names = F)
sel <- sel[sel$Qvalue < 0.05,]
qvalue <- qvalue[unique(sel$X1),unique(sel$X2)]
cor <- cor[unique(sel$X1),unique(sel$X2)]
qvalue <- qvalue[,-grep('Cactin',colnames(qvalue))]
cor <- cor[,-grep('Cactin',colnames(cor))]

annotation_col <- data.frame(pathway=c(rep('NAD+ Metabolism',12),rep('Fatty Acid Metabolism',8)))
rownames(annotation_col) <- c('Cd73','Nadk','Nadsyn1','Nampt','Naprt','Nmnat3','Nnt','Nrk1',
                              'Paps','Pnp','Qprt','Sirts','Acadl','Acsl4','Cact','Cd36','Cpt1a','Cpt2','Crat','Fatp2')
cor <- cor[,c(1,3,4,5,6,7,2,8,9,10,11,12)]
qvalue <- qvalue[,c(1,3,4,5,6,7,2,8,9,10,11,12)]
pathway <- jet(2)
names(pathway) <- unique(annotation_col[,1])
annotation_color <- list(pathway=pathway)
limits<-floor(max(max(cor,na.rm = T),1*min(cor,na.rm = T))*100000)/100000
bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
pheatmap(cor,scale = "none",color = c,show_colnames=T,cluster_cols=F,cluster_rows = F,display_numbers = ifelse(qvalue<0.05,"*",""),annotation_col = annotation_col,#annotation_row = annotation,
         breaks=1*bk,cellwidth = 20, cellheight = 10, fontsize_row  = 10, fontsize_col = 10, border_color = "white",show_rownames = T,fontsize_number = 18,
         filename = paste0(k," 相关性heatmap.png"),annotation_colors = annotation_color)
pheatmap(cor,scale = "none",color = c,show_colnames=T,cluster_cols=F,cluster_rows = F,display_numbers = ifelse(qvalue<0.05,"*",""),annotation_col = annotation_col,#annotation_row = annotation,
         breaks=1*bk,cellwidth = 20, cellheight = 10, fontsize_row  = 10, fontsize_col = 10, border_color = "white",show_rownames = T,fontsize_number = 18,
         filename = paste0(k," 相关性heatmap.pdf"),annotation_colors = annotation_color)
##############
tfgene <- read.xlsx('06.TF/线粒体nad_fatty_Clock.xlsx',sheet = 1,colNames = T,rowNames = F,check.names = F)
tflist <- unique(c(tfgene$TF,tfgene$Gene))
data1 <- as.data.frame(t(allrhy[tflist,grep('DM',colnames(allrhy))]))
data2 <- as.data.frame(t(allrhy[tflist,grep('DM',colnames(allrhy))]))

for(i in 1:nrow(tfgene)){
  tfgene[i,3] <- cor[tfgene[i,1],tfgene[i,2]]
  tfgene[i,4] <- pvalue[tfgene[i,1],tfgene[i,2]]
}
write.xlsx(tfgene,'tfgene_network.xlsx')

tfgene$TF[match(unique(tfgene$Gene),tfgene$TF)]
dmrhy <- read.csv('01.节律与差异节律/dm_rhythm.csv',row.names = 1)
rownames(ctrl)[match(colnames(cor),rownames(dmrhy[dmrhy$pVal<0.05,]))]
