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
pdf('mfuzz.pdf',width = 8,height = 6)
mfuzz.plot2(df3F, cl=cl,mfrow=c(2,3),centre=TRUE,x11=F,centre.lwd=0.2)
dev.off()
#批量导出每个聚类所包含的基因
dir.create(path="03.mfuzz",recursive = TRUE)
for(i in 1:6){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("03.mfuzz","/mfuzz_",i,".csv"))
}
