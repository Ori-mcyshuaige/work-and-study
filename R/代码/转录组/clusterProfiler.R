library(org.Hs.eg.db) #人类注释数据库
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(openxlsx)
library(org.Mm.eg.db)

genelist <- read.table('ML/00.data/merge.txt',header = T,row.names = 1)
genelist <- rownames(genelist)

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
go@result <- go@result[order(go@result$FoldEnrichment,decreasing = T),]
write.csv(go,file="go.csv")
bar <- barplot(go, split="ONTOLOGY",showCategory=10,x='FoldEnrichment')+ facet_grid(ONTOLOGY~.,scale="free")
ggsave('go_barplot.pdf',bar,dpi = 600,height = 12,width = 10)


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
kegg@result <- kegg@result[order(kegg@result$FoldEnrichment,decreasing = T),]
write.csv(kegg,file = "kegg.csv")
kegg_bar <- barplot(kegg,showCategory=20,drop=T) #柱状图
ggsave('kegg_barplot.pdf',kegg_bar,dpi = 600,height = 12,width = 10)
