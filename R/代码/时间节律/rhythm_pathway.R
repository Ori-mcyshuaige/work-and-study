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
options(clusterProfiler.download.method = "wininet")
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
#小鼠所有的基因集
m_df = msigdbr(species = "Mus musculus")
##2.3 查看基因集类别：
a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
View(a)
##2.5 检索鼠类C2 (curated) CGP (chemical and genetic perturbations)基因集：
m_df = msigdbr(species = "Mus musculus", category = "C2" )
##KEGG   
m_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
#C5   GO:BP CC MF  HPO
m_df = msigdbr(species = "Mus musculus", category = "C5")
##2.4 检索鼠类的hallmark 基因集：
m_df = msigdbr(species = "Mus musculus", category = "H")
## 免疫相关通路，C7
m_df = msigdbr(species = "Mus musculus", category = "C7")

m_df1 =m_df[,c("gs_name","gs_url")] 
m_df1 = unique(m_df1)
m_df1$gs_symbol = NA
for(i in m_df1$gs_name){
  m_df1[m_df1$gs_name==i,'gs_symbol']=paste0(unlist(m_df[m_df$gs_name==i,"gene_symbol"]),collapse = "\t")
}
write.table(m_df1,'mh.cp2.kegg.gmt',sep = "\t",row.names = F,col.names = F,quote = F)

ctrl <- read.csv('01.节律与差异节律/ctrl_rhygene.csv',row.names = 1)
dm <- read.csv('01.节律与差异节律/dm_rhygene.csv',row.names = 1)
allrhy <- cbind(ctrl,dm[rownames(ctrl),])
allrhy <- read.csv('01.节律与差异节律/difrhygene_log2cpm+1.csv',row.names = 1)

geneSets=getGmt('./mh.cp2.kegg.gmt', geneIdType=SymbolIdentifier())
ssgseaScore=gsva(as.matrix(allrhy), geneSets, method='gsva')

#对GSVA的打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
write.csv(ssgseaScore,'gsva_res.csv')
Type1=sapply(colnames(ssgseaScore),function(x){
  strsplit(x,'_')[[1]][1]
})
Type2=sapply(colnames(ssgseaScore),function(x){
  strsplit(x,'_')[[1]][2]
})
#通路差异分析
outTab=data.frame()
for(i in row.names(ssgseaScore)){
  test=t.test(ssgseaScore[i,] ~ Type1)
  pvalue=test$p.value
  t=test$statistic
  log2fc=log2(mean(ssgseaScore[i,grep('DM',colnames(ssgseaScore))])/mean(ssgseaScore[i,grep('Ctrl',colnames(ssgseaScore))]))
  if(pvalue<0.05){
    Sig=ifelse(pvalue>0.05, "Not", ifelse(log2fc>0,"Up","Down"))
    outTab=rbind(outTab, cbind(Pathway=i, t,log2fc, pvalue, Sig))
  }
}
write.csv(outTab,'gsva_stat_p.csv')
###plot heatmap
# data <- ssgseaScore[outTab$Pathway,]
data <- ssgseaScore[rownames(keepgene),]
for(ii in 1:nrow(data)){
  result<-sapply(as.numeric(data[ii,]),function(i){(i-mean(as.numeric(data[ii,])))/sd(as.numeric(data[ii,]))})
  data[ii,]<-result
}

limits<-floor(max(max(data,na.rm = T),1*min(data,na.rm = T))*100000)/100000
bk <- c(seq(-1*limits,-0.001,by=0.001),seq(0,limits,by=0.001))
c<-c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

Cluster <- data.frame(Group=Type1,Time=Type2)
rownames(Cluster) <- colnames(data)
ann_colors=list()
Type1c=jet(2)
Type2c=jet(6)
names(Type1c)=levels(factor(Cluster$Group))
names(Type2c)=levels(factor(Cluster$Time,levels = unique(Cluster$Time)))
ann_colors[["Group"]]=Type1c
ann_colors[["Time"]]=Type2c
##pheatmap函数绘制热图
pheatmap(data,scale = "none",color = c,show_colnames=F,show_rownames = T,cluster_cols=F,cluster_rows = T,annotation_col = Cluster,#clustering_method = "complete",#annotation_row = annotation_row,
         breaks=1*bk,cellwidth = 6, cellheight = 10, fontsize_row  = 10, fontsize_col = 20, border_color = "white",annotation_colors = ann_colors,
         filename = "heatmap.pdf")

#绘制柱状图
termNum=5      #展示通路的数目
outTab=outTab[order(outTab$t),]
outTab=outTab[c(1:termNum,(nrow(outTab)-termNum+1):nrow(outTab)),]
pdf(file="./barplot.pdf", width=9, height=6)
outTab$log2fc=as.numeric(outTab$log2fc)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="log2fc", fill = "Sig", color = "white",
              palette=c("#0066FF","#FF9900"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title="",
              xlab="Term", ylab="log2fc(GSVA score, DM/Ctrl)",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()


#######################
raw <- read.csv('01.节律与差异节律/difrhygene_log2cpm+1.csv',row.names = 1)
raw <- as.data.frame(t(raw))
raw$group <- sapply(rownames(raw),function(x){paste0(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2])})
df2<-aggregate(raw[,-ncol(raw)],by=list(raw$group),mean,na.rm= TRUE)
rownames(df2) <- df2$Group.1
df2 <- df2[,-1]
for(i in colnames(df2)){
  df2[,i] <- as.numeric(df2[,i])
}
write.csv(df2,'mean_dfrhygene.csv')

out <- data.frame()
for(i in colnames(df2)){
  Ctrl <- rownames(df2)[which.max(df2[grep('Ctrl',rownames(df2)),i])]
  DM <- rownames(df2)[6+which.max(df2[grep('DM',rownames(df2)),i])]
  out <- rbind(out,cbind(gene=i,Ctrl,DM))
}

geneset <- list()
for(i in rownames(df2)){
  geneset[[i]] <- unlist(out$gene[out$Ctrl==i|out$DM==i])
}
# for(i in colnames(data)){
#   geneset[[i]] <- unlist(data[,i][!is.na(data[,i])])
# }


GO_BP <- data.frame()
GO_CC <- data.frame()
GO_MF <- data.frame()
KEGG <- data.frame()
Reactome <- data.frame()
for(i in names(geneset)){
  genelist <- geneset[[i]]
  
  hg <- bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),
             OrgDb="org.Mm.eg.db")## 鼠 org.Mm.eg.db   人org.Hs.eg.db
  head(hg)
  go <- enrichGO(hg$ENTREZID,OrgDb = org.Mm.eg.db, ont='ALL',
                 pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID')#进行GO富集，确定P值与Q值得卡值并使用BH方法对值进行调整。
  go@result$FoldEnrichment <- apply(go@result,1,function(x){
    GeneRatio=eval(parse(text=x["GeneRatio"]))
    BgRatio=eval(parse(text=x["BgRatio"]))
    enrichment_fold=round(GeneRatio/BgRatio,2)
    enrichment_fold
  })
  go@result$Time <- strsplit(i,'_')[[1]][2]
  go@result$Group <- strsplit(i,'_')[[1]][1]
  go@result <- go@result[order(go@result$FoldEnrichment,decreasing = T),]
  write.csv(go,paste0(i,'_go.csv'))
  kegg <- enrichKEGG(hg$ENTREZID, organism = 'mmu', keyType = 'kegg',
                     pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 3,
                     maxGSSize = 500, qvalueCutoff = 0.2,
                     use_internal_data = FALSE)#进行KEGG富集 鼠mmu  人 hsa
  kegg <- setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  kegg@result$FoldEnrichment <- apply(kegg@result,1,function(x){
    GeneRatio=eval(parse(text=x["GeneRatio"]))
    BgRatio=eval(parse(text=x["BgRatio"]))
    enrichment_fold=round(GeneRatio/BgRatio,2)
    enrichment_fold
  })
  kegg@result$Time <- strsplit(i,'_')[[1]][2]
  kegg@result$Group <- strsplit(i,'_')[[1]][1]
  kegg@result <- kegg@result[order(kegg@result$FoldEnrichment,decreasing = T),]
  write.csv(kegg,paste0(i,'_kegg.csv'))
  react = ReactomePA::enrichPathway(gene = hg$ENTREZID,organism = 'mouse')
  react <- setReadable(react, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  react@result$FoldEnrichment <- apply(react@result,1,function(x){
    GeneRatio=eval(parse(text=x["GeneRatio"]))
    BgRatio=eval(parse(text=x["BgRatio"]))
    enrichment_fold=round(GeneRatio/BgRatio,2)
    enrichment_fold
  })
  react@result$Time <- strsplit(i,'_')[[1]][2]
  react@result$Group <- strsplit(i,'_')[[1]][1]
  react@result <- react@result[order(react@result$FoldEnrichment,decreasing = T),]
  write.csv(react,paste0(i,'_react.csv'))
  
  bp <- go@result[go@result$ONTOLOGY=='BP'&go@result$p.adjust<0.05,c('Description','Time','p.adjust','FoldEnrichment','Group','geneID')]
  cc <- go@result[go@result$ONTOLOGY=='CC'&go@result$p.adjust<0.05,c('Description','Time','p.adjust','FoldEnrichment','Group','geneID')]
  mf <- go@result[go@result$ONTOLOGY=='MF'&go@result$p.adjust<0.05,c('Description','Time','p.adjust','FoldEnrichment','Group','geneID')]
  kegg <- kegg@result[kegg@result$p.adjust<0.05,c('Description','Time','p.adjust','FoldEnrichment','Group','geneID')]
  react <- react@result[react@result$p.adjust<0.05,c('Description','Time','p.adjust','FoldEnrichment','Group','geneID')]
  if(nrow(bp)>0){
    GO_BP <- rbind(GO_BP,bp[1:ifelse(nrow(bp)<10,nrow(bp),10),])
  }
  if(nrow(cc)>0){
    GO_CC <- rbind(GO_CC,cc[1:ifelse(nrow(cc)<10,nrow(cc),10),])
  }
  if(nrow(mf)>0){
    GO_MF <- rbind(GO_MF,mf[1:ifelse(nrow(mf)<10,nrow(mf),10),])
  }
  if(nrow(kegg)>0){
    KEGG <- rbind(KEGG,kegg[1:ifelse(nrow(kegg)<10,nrow(kegg),10),])
  }
  if(nrow(react)>0){
    Reactome <- rbind(Reactome,react[1:ifelse(nrow(react)<10,nrow(react),10),])
  }
}


GO_BP <- GO_BP[c(grep('ZT0',GO_BP$Time),grep('ZT4',GO_BP$Time),
               grep('ZT8',GO_BP$Time),grep('ZT12',GO_BP$Time),
               grep('ZT16',GO_BP$Time),grep('ZT20',GO_BP$Time)),]
GO_CC <- GO_CC[c(grep('ZT0',GO_CC$Time),grep('ZT4',GO_CC$Time),
               grep('ZT8',GO_CC$Time),grep('ZT12',GO_CC$Time),
               grep('ZT16',GO_CC$Time),grep('ZT20',GO_CC$Time)),]
GO_MF <- GO_MF[c(grep('ZT0',GO_MF$Time),grep('ZT4',GO_MF$Time),
               grep('ZT8',GO_MF$Time),grep('ZT12',GO_MF$Time),
               grep('ZT16',GO_MF$Time),grep('ZT20',GO_MF$Time)),]
KEGG <- KEGG[c(grep('ZT0',KEGG$Time),grep('ZT4',KEGG$Time),
               grep('ZT8',KEGG$Time),grep('ZT12',KEGG$Time),
               grep('ZT16',KEGG$Time),grep('ZT20',KEGG$Time)),]
Reactome <- Reactome[c(grep('ZT0',Reactome$Time),grep('ZT4',Reactome$Time),
               grep('ZT8',Reactome$Time),grep('ZT12',Reactome$Time),
               grep('ZT16',Reactome$Time),grep('ZT20',Reactome$Time)),]
GO_BP$Time <- factor(GO_BP$Time,levels = unique(GO_BP$Time))
GO_BP$Description <- factor(GO_BP$Description,levels = unique(GO_BP$Description))
GO_CC$Time <- factor(GO_CC$Time,levels = unique(GO_CC$Time))
GO_CC$Description <- factor(GO_CC$Description,levels = unique(GO_CC$Description))
GO_MF$Time <- factor(GO_MF$Time,levels = unique(GO_MF$Time))
GO_MF$Description <- factor(GO_MF$Description,levels = unique(GO_MF$Description))
KEGG$Time <- factor(KEGG$Time,levels = unique(KEGG$Time))
KEGG$Description <- str_split_i(KEGG$Description,' - ',1)
KEGG$Description <- factor(KEGG$Description,levels = unique(KEGG$Description))
Reactome$Time <- factor(Reactome$Time,levels = unique(Reactome$Time))
Reactome$Description <- str_split_i(Reactome$Description,' - ',1)
Reactome$Description <- factor(Reactome$Description,levels = unique(Reactome$Description))


pdata <- list(GO_BP=GO_BP,GO_CC=GO_CC,GO_MF=GO_MF,KEGG=KEGG,Reactome=Reactome)
for(k in names(pdata)){
  data <- pdata[[k]]
  p <- ggplot(data,aes(x=Time,y=Description,fill=p.adjust,
                       color=p.adjust,size=FoldEnrichment))+
    geom_point(shape=16)+
    # scale_x_discrete(limits = c('ZT0','ZT4','ZT8','ZT12','ZT16','ZT20'))+
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
      name = 'FoldEnrichment', guide = guide_legend(order = 2), range = c(3, 9),
      breaks = c(min(abs(data$FoldEnrichment)), max(abs(data$FoldEnrichment))),
      labels = c(round(min(data$FoldEnrichment), 1), round(max(data$FoldEnrichment), 1))
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
    facet_grid(~data$Group,scales= "free",space= "free")+
    labs(x = NULL,y = k,color=NULL) 
  
  ggsave(
    paste0(k,"_",'all.pdf'), p,
    width = 12, height = 10, units = 'in', dpi = 600
  )
  ggsave(
    paste0(k,"_",'all.png'), p,
    width = 12, height = 10, units = 'in'
  )
}

##########################
ampdata <- data.frame(matrix(nrow = 12,ncol = 3))
colnames(ampdata) <- c('Group','Time','number')
ampdata$Time <- str_remove(ampdata$Time,'ZT')
for(i in 1:12){
  ampdata[i,1] <- strsplit(names(geneset)[i],'_')[[1]][1]
  ampdata[i,2] <- strsplit(names(geneset)[i],'_')[[1]][2]
  ampdata[i,3] <- length(geneset[[i]])
}
ampdata$percentage <- ampdata$number/502
p<-ggplot(data = ampdata,aes(x = Time, y = percentage,color=Group,fill=Group))+
  geom_bar(stat = "identity",position = "dodge",width=.8)+  ###position = "dodge" 列  “stack”堆叠   "fill" 百分比
  scale_x_discrete(limits = c('0','4','8','12','16','20'))+
  scale_color_manual(values = c("DM"="red","Ctrl"="black"))+
  # geom_text(aes(x=label,y=Tissue,label=Number),size=7,color="black",position='identity')+
  scale_fill_manual(values = c("DM"="red","Ctrl"="black"))+
  theme_bw()+
  theme(
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",size = 1),
    legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_blank(),
    # legend.position = "none",
    # legend.key.size = unit(0.6,"cm"),
    axis.title = element_text(size = 14,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text.x = element_text(size = 14,color = "black",family = "sans",face = "plain",angle = 0, hjust = 0,vjust = 0),
    axis.text.y = element_text(size = 14,color = "black",family = "sans",face = "plain", hjust = 0,vjust = 0)
  )+
  # scale_y_continuous(limits = c(0,58),expand = c(0,0))+
  scale_y_continuous(labels = scales::percent_format())+
  coord_polar()+
  labs(x = "Time(hours)", y = "Percentage")
ggsave("percentage.png",p,dpi = 600,width = 6,height = 6)
ggsave("percentage.pdf",p,width = 6,height = 6,dpi = 600)


clkgdata <- data.frame(matrix(nrow = 11,ncol = 2))
clkgs <- c('Ciart','Clock','Cry1','Cry2','Dbp','Npas2','Nr1d1','Nr1d2','Per1','Per2','Per3')
colnames(clkgdata) <- c('Ctrl','DM')
rownames(clkgdata) <- clkgs
for (i in clkgs){
  clkgdata[i,'Ctrl'] <- strsplit(names(geneset[1:6])[grep(i,geneset[1:6])],'_')[[1]][2]
  clkgdata[i,'DM'] <- strsplit(names(geneset[7:12])[grep(i,geneset[7:12])],'_')[[1]][2]
}
# clkgdata$Ctrl <- as.numeric(str_remove(clkgdata$Ctrl,'ZT'))
# clkgdata$DM <- as.numeric(str_remove(clkgdata$DM,'ZT'))
clkgdata <- as.data.frame(as.table(as.matrix(clkgdata)))
colnames(clkgdata) <- c('clkgs','group','time')
clkgdata$time <- str_remove(clkgdata$time,'ZT')
clkgdata$time <- factor(clkgdata$time,levels = c('0','4','8','12','16','20'))
Ctrl <- read.csv('01.节律与差异节律/ctrl_rhythm.csv',row.names = 1)
DM <- read.csv('01.节律与差异节律/dm_rhythm.csv',row.names = 1)
clkgdata$p <- c(Ctrl[clkgs,]$pVal,DM[clkgs,]$pVal)
clkgdata$psig <- ifelse(clkgdata$p<0.001,'***',ifelse(clkgdata$p<0.01,'**',ifelse(clkgdata$p<0.05,'*',"")))
p<-ggplot(clkgdata,aes(clkgs,group))+
  geom_tile(color="white",aes(fill=time),size=1)+
  geom_text(aes(x = clkgs,y = group,label=psig),size=10,color="black")+
  scale_fill_manual(
    name = NULL, values = c(
      '0'="#30123B", '4'="#3D9BFC", '8'="#47F682", '12'="#E0DB37",'16'="#EE5A11", '20'="#7A0403"
    ), guide = guide_legend(order = 1, override.aes = list(size = 3),color='black')
  )+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(family="sans"),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,size = 15,color='black'),
        axis.text.y = element_text(size = 15,color='black'),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill=NA,color="white",
                                    size=2, linetype="solid"))
ggsave(
  'clkgene_amp.jpg', p,
  width = 9, height = 2, units = 'in', dpi = 600
)
ggsave(
  'clkgene_amp.pdf', p,
  width = 9, height = 2, units = 'in'
)

