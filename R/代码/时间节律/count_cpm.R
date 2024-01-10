library(edgeR)

counts <- read.table('00.data/result/gene_counts/gene_counts.xls',row.names = 1,header = T)
cpm <- edgeR::cpm(counts)

out <- cpm[sapply(rownames(cpm),function(x){all(mean(unlist(cpm[x,grep('ZT0',colnames(cpm))]))>0.5,
                                     mean(unlist(cpm[x,grep('ZT4',colnames(cpm))]))>0.5,
                                     mean(unlist(cpm[x,grep('ZT8',colnames(cpm))]))>0.5,
                                     mean(unlist(cpm[x,grep('ZT12',colnames(cpm))]))>0.5,
                                     mean(unlist(cpm[x,grep('ZT16',colnames(cpm))]))>0.5,
                                     mean(unlist(cpm[x,grep('ZT20',colnames(cpm))]))>0.5
                                     )}),]
out <- log2(out+1)
write.csv(out,'keeplogcpm.csv')

cpm <- cpm[edgeR::filterByExpr(cpm),]
