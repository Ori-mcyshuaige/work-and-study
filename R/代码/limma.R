library(limma)

pred='group'
con='control'
disease='MDD'
data <- t(read.csv('./01.差异分析/GSE98793.csv',header = T,row.names = 1))


design <- as.data.frame(model.matrix(~0+factor(data[,pred])))
rownames(design) <- rownames(data)
colnames(design) <- unique(data[,pred])
data <- as.data.frame(t(data[,-match(pred,colnames(data))]))
# colnames(data) <- rownames(data[,-match(pred,colnames(data))])
# rownames(data) <- colnames(data[,-match(pred,colnames(data))])
for(i in colnames(data)){
  data[,i] <- as.numeric(data[,i])
}
fit <- lmFit(data,design)
cont.matrix <- makeContrasts(paste0(disease,'-',con),levels=design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
genestat <- topTable(fit2,adjust = 'fdr',number=200000)
