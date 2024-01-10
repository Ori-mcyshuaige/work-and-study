library(tidyverse)
library(tidyr)
library(dplyr)
library(reprex)
library(stringr)

library(rtracklayer)
gtf <- rtracklayer::import('Homo_sapiens.GRCh38.109.gtf')
gtf<-as.data.frame(gtf)
##获取编码蛋白gene
protein_coding.gtf<-gtf %>% 
dplyr::select(c("gene_id","gene_name","gene_biotype","gene_biotype")) %>% 
arrange(gene_id) %>% 
distinct(gene_name,.keep_all = T) %>% 
filter(gene_biotype == "protein_coding")

##同理可以获取lnRNA、miRNA
lncRNA.gtf<-gtf %>% 
dplyr::select(c("gene_id","gene_name","gene_biotype","gene_biotype")) %>% 
arrange(gene_id) %>% 
distinct(gene_name,.keep_all = T) %>% 
filter(gene_biotype == "lncRNA")

miRNA.gtf<-gtf %>% 
dplyr::select(c("gene_id","gene_name","gene_biotype","gene_biotype")) %>% 
arrange(gene_id) %>% 
distinct(gene_name,.keep_all = T) %>% 
filter(gene_biotype == "miRNA")

write.csv(lncRNA.gtf,'lncRNA_db.csv',row.names = F)
write.csv(miRNA.gtf,'miRNA_db.csv',row.names = F)
write.csv(protein_coding.gtf,'mRNA_db.csv',row.names = F)
