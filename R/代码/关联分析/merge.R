
setwd("C:\\Users\\Administrator\\Desktop\\关联分析\\数据整理")
library("openxlsx")
aa<-read.xlsx("Metagenomics.xlsx")
bb<-read.xlsx("Metagenomicsd.xlsx")

cc <- merge(aa,bb,by='OTU')
write.xlsx(cc,paste("meta.xlsx"))