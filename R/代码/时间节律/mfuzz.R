library(Mfuzz)

raw <- read.csv('01.节律与差异节律/keeplogcpm.csv',row.names = 1)
dmrhy <- read.csv('01.节律与差异节律/dm_rhythm.csv',row.names = 1)
dm <- raw[rownames(dmrhy)[dmrhy$pVal < 0.05],]
dm <- dm[,grep('DM',colnames(dm))]
dm <- dm[,c(grep('ZT0',colnames(dm)),grep('ZT4',colnames(dm)),grep('ZT8',colnames(dm)),grep('ZT12',colnames(dm)),
            grep('ZT16',colnames(dm)),grep('ZT20',colnames(dm)))]
dm <- as.data.frame(t(dm))
dm$group <- sapply(rownames(dm),function(x){strsplit(x,'_')[[1]][2]})
df2<-aggregate(dm[,-ncol(dm)],by=list(dm$group),mean,na.rm= TRUE)
rownames(df2) <- df2[,1]
df3<-t(df2[,-1])
df3Ex<- ExpressionSet(assayData = df3)

#排除了超过25%的测量缺失的基因
df3F <- filter.NA(df3Ex,thres = 0.25)
## 过滤标准差为0的基因
df3F <- filter.std(df3F,min.std=0)
#用相应基因的平均值表达值替换剩余的缺失值
df3F <- fill.NA(df3F,mode = 'mean')
#标准化
df3F <- standardise(df3F)


## 聚类个数
c <- 6
## 计算最佳的m值
m <- mestimate(df3F)
## 聚类
cl <- mfuzz(df3F, c = c, m = m)


## 查看每类基因数目
cl$size
## 查看每类基因ID
cl$cluster[cl$cluster == 1]
## 输出基因ID
# write.table(cl$cluster,"output.txt",quote=F,row.names=T,col.names=F,sep="\t")
## 绘制折线图
# mfuzz.plot(df3F,cl,mfrow=c(2,3),new.window= FALSE)
mfuzz.plot2(df3F, cl=cl,mfrow=c(2,3),centre=TRUE,x11=F,centre.lwd=0.2)
#批量导出每个聚类所包含的基因
dir.create(path="mfuzz",recursive = TRUE)

for(i in 1:6){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("03.mfuzz","/mfuzz_",i,".csv"))
}

geneset <- list()
for(i in 1:6){
  cls <- paste0('cluster',i)
  geneset[[cls]] <-names(cl$cluster[cl$cluster == i])
}


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
  p <- ggplot(data,aes(x=Group,y=Description,fill=p.adjust,
                       color=p.adjust,size=FoldEnrichment))+
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
    # facet_grid(~data$Group,scales= "free",space= "free")+
    labs(x = NULL,y = k,color=NULL) 
  
  ggsave(
    paste0(k,"_",'all.pdf'), p,
    width = 15, height = 10, units = 'in', dpi = 600
  )
  ggsave(
    paste0(k,"_",'all.png'), p,
    width = 15, height = 10, units = 'in'
  )
}
