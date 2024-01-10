library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
library(rio)
library(GSEABase)

data <- read.csv('01.节律与差异节律/keeplogcpm.csv',row.names = 1)
data <- data[,c(grep('ZT0',colnames(data)),grep('ZT4',colnames(data)),grep('ZT8',colnames(data)),grep('ZT12',colnames(data)),
                      grep('ZT16',colnames(data)),grep('ZT20',colnames(data)))]
data <- data[,c(grep('Ctrl',colnames(data)),grep('DM',colnames(data)))]

nadmet <- c('Nadk','Nadsyn1','Nampt','Naprt','Nmnat3','Nnt','Pnp')
fattymet <- c('Acadl','Acsl4','Cd36','Cpt1a','Cpt2','Crat')
# match(fattymet,rownames(data))
geneset <- list(NADMET=nadmet,FattyMET=fattymet)
ssgseaScore=gsva(as.matrix(data), geneset, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
write.csv(ssgseaScore,file="ssGSEA.metscore.csv",row.names = T)
