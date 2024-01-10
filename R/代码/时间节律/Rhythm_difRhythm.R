if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rain")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("npsm")

library(devtools)
library(rain)
library(DODR)

raw <- read.csv('./01.节律与差异节律/keeplogcpm.csv',row.names = 1)
control <- raw[,grep('Ctrl',colnames(raw))]
control <- control[,c(grep('ZT0',colnames(control)),grep('ZT4',colnames(control)),grep('ZT8',colnames(control)),grep('ZT12',colnames(control)),
                      grep('ZT16',colnames(control)),grep('ZT20',colnames(control)))]
dm <- raw[,grep('DM',colnames(raw))]
dm <- dm[,c(grep('ZT0',colnames(dm)),grep('ZT4',colnames(dm)),grep('ZT8',colnames(dm)),grep('ZT12',colnames(dm)),
                      grep('ZT16',colnames(dm)),grep('ZT20',colnames(dm)))]


timedm <- sapply(colnames(dm),function(x){paste0(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2])})
timecontrol <- sapply(colnames(control),function(x){paste0(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2])})

ctrl_rs<- rain(t(control),
               deltat = 4, #时间点间隔
               period = 24,#周期
               nr.series = 4, #每一时间点重复次数
               measure.sequence = c(4,4,3,4,4,4), #单独指定每一时间点重复次数，允许某一时间点没有测量，会覆盖前面nr.series参数
               method = "independent",#检测方法，"independent"代表多个周期被解释为一个周期的重复 ；"logitudinal"代表进行的是纵向数据（重复测量，不规则采样，不相等的时间间隔。）
               peak.border = c(0.3, 0.7),#定义峰型，一般默认
               adjp.method = "ABH", 
               verbose=TRUE)
write.csv(ctrl_rs,'./ctrl_rhythm.csv')

dm_rs<- rain(t(dm),
               deltat = 4, #时间点间隔
               period = 24,#周期
               nr.series = 4, #每一时间点重复次数
               measure.sequence = NULL, #单独指定每一时间点重复次数，允许某一时间点没有测量，会覆盖前面nr.series参数
               method = "independent",#检测方法，"independent"代表多个周期被解释为一个周期的重复 ；"logitudinal"代表进行的是纵向数据（重复测量，不规则采样，不相等的时间间隔。）
               peak.border = c(0.3, 0.7),#定义峰型，一般默认
               adjp.method = "ABH", 
               verbose=TRUE)
write.csv(dm_rs,'./dm_rhythm.csv')

samegene <- intersect(rownames(ctrl_rs)[ctrl_rs$pVal < 0.05],rownames(dm_rs)[dm_rs$pVal < 0.05])
control <- control[samegene,]
dm <- dm[samegene,]
write.csv(control,'control_rhygene.csv')
write.csv(dm,'dm_rhygene.csv')


########difRhythm,dodr包：
# dodr(val1, #行为时间点，列为基因名的矩阵1
#      val2, #行为时间点，列为基因名的矩阵2
#      times1, #时间点序列1
#      times2 = times1, #时间点序列2，默认等于1
#      norm = TRUE, #是否在分析之前对时间序列进行归一化(按平均值分割)
#      period = 24,#周期
#      method = "robust", #方法，默认为“robust”，可选有6种
#      verbose = options("verbose")[[1]])
#示例
# n=50
# testTimes1 <- 0:15*3
# testTimes2 <- testTimes1
# tp <- length(testTimes1)
# per1 <- 24
# amp1 <- 0.3
# ph1 <- 5
# sd1 <- 0.1
# 
# per2 <- per1
# amp2 <- amp1
# ph2 <- ph1+4
# sd2 <- sd1
# 
# #creating artificial oscillation sets
# v1 <- 1 + amp1 * cos((testTimes1 - ph1)/per1*2*pi)
# noise1 <- rnorm(length(testTimes1)*n, 0, sd1)
# val1 <- matrix(v1 + noise1, ncol=n)
# 
# v2 <- 1 + amp2 * cos((testTimes2 - ph2)/per2*2*pi)
# noise2 <- rnorm(length(testTimes2)*n, 0, sd2)
# val2 <- matrix(v2 + noise2, ncol=n)
# 
# # run DODR
# dodr <- dodr(val1, val2, testTimes1, testTimes2, 24, method = 'all')
# dodr$p.value.table
val1 <- t(control)
val2 <- t(dm)
times1 <- rep(0:5*4,each=4)
times2 <- times1
times1 <- times1[-9]
difrhy <- dodr(val1, #行为时间点，列为基因名的矩阵1
     val2, #行为时间点，列为基因名的矩阵2
     times1 = times1, #时间点序列1
     times2 = times2, #时间点序列2，默认等于1
     norm = TRUE, #是否在分析之前对时间序列进行归一化(按平均值分割)
     period = 24,#周期
     method = "all", #方法，默认为“robust”，可选有6种
     verbose = options("verbose")[[1]])
difrhy_rs <- difrhy$p.value.table
rownames(difrhy_rs) <- rownames(control)
write.csv(difrhy_rs,'difrhy_rs.csv')
keepgene <- cbind(control,dm)[rownames(difrhy_rs)[difrhy_rs$meta.p.val<0.05],]
write.csv(keepgene,'difrhygene_log2cpm+1.csv')




