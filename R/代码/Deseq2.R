library(DESeq2)
library(dplyr)
mydata3 <- read.table('GSE236761_raw_counts.txt',header = T,row.names = 1)
coldata1 <- data.frame(condition=factor(c(rep('ASD',12),rep('Ctrl',8)),
                                        levels = c('ASD','Ctrl')))

dds1 <- DESeqDataSetFromMatrix(countData = as.data.frame(mydata3), colData = coldata1, design= ~condition)
dds1.1 <- DESeq(dds1)
v = vst(dds1.1, blind=FALSE)
vsted = assay(v)#######该数据可作为CBNplot的输入数据
res <- results(dds1.1, contrast = c('condition','ASD','Ctrl'))
res1 <- data.frame(res,stringsAsFactors = F,check.names = F)
res1 <- na.omit(res1)
data1 <- res1 %>% 
  mutate(change=as.factor(ifelse(res1$pvalue<0.05 & abs(res1$log2FoldChange)>log2(2),
                                 ifelse(res1$log2FoldChange>0,'up','down'),
                                 'no change')))
sig <- data1[data1$change %in% c('up','down'),]
vsted <- vsted[c(rownames(sig),'TSC1'),]
write.csv(sig,'dif_gene.csv')
