# devtools::install_github("IOBR/IOBR",ref="master")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

rm(list = ls())
library(dplyr)
library(IOBR)
library(UCSCXenaTools)
library(limma)
help("count2tpm")

# eset_prad<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Prostate Cancer (PRAD)") %>% 
#   XenaFilter(filterDatasets    = "TCGA-PRAD.htseq_counts.tsv") %>% 
#   XenaQuery() %>%
#   XenaDownload() %>% 
#   XenaPrepare()
# eset_prad$Ensembl_ID <- substring(eset_prad$Ensembl_ID,1,15)
# eset_prad <- column_to_rownames(eset_prad,var = "Ensembl_ID")
# eset_prad<-(2^eset_prad)-1




mycounts <- as.matrix(read.table('./ML/00.data/RNA_matrix.txt',header = T))
mycounts <- mycounts[mycounts[,"Geneid"]!='',]


# mycounts<-read.csv("merge.csv")
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)#最后一列Length是基因长度
mycounts <- avereps(mycounts)
mycounts <- as.data.frame(mycounts)
data <- as.data.frame(lapply(mycounts,as.numeric))
rownames(data) <- rownames(mycounts)
#Counts转TPM

tpm <- count2tpm(countMat = data,idType = "SYMBOL",source = "local")
write.csv(tpm,'ourdata_tpm.csv')

Counts2TPM <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

#Counts转FPKM

Counts2FPKM <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

#fpkm转tpm
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}