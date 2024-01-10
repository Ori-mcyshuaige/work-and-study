rm(list = ls())
options(stringsAsFactors = F)

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("GEOquery")
BiocManager::install("AnnoProbe")

library(TCGAbiolinks)
library(scRNAseq)
library(data.table)
library(limma)
library(dplyr)
library(DT)
library(survival)
library(survminer)
query <- GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")
GDCdownload(query = query)
expData <- GDCprepare(query = query,
                      save = T,
                      save.filename = "BLCA_exp.rda")
tpm_data <- assay(expData,i = "tpm_unstrand")
row_file <- data.table::fread("./GDCdata/TCGA-BLCA/Transcriptome_Profiling/Gene_Expression_Quantification/0aacd6c6-c2cb-4304-8d2d-ff63e15629c2/cb8fa4e8-736f-4e98-98a7-5d6ea4dcb800.rna_seq.augmented_star_gene_counts.tsv",
                              data.table = F)
row_file <- row_file[-c(1:4),]
rownames(row_file) <- row_file[,1]
same <- intersect(row.names(tpm_data),row.names(row_file))
length(same)
BLCA_tpm <- cbind(row_file[same,],tpm_data[same,])
BLCA_tpm <- BLCA_tpm[,-c(1,3:9)]
dim(BLCA_tpm)
#去重
rt <- as.matrix(BLCA_tpm)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <-list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
BLCA_exp <- data[rowMeans(data)>0,]
dim(BLCA_exp)
Out=rbind(id=colnames(BLCA_exp), BLCA_exp)
write.table(Out, file="./survival/BLCA_exp.txt", sep="\t", quote=F, col.names=F)
####下载临床数据
query <- GDCquery(
  project = "TCGA-BLCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)
GDCdownload(query)
cli <- GDCprepare_clinic(query,'follow_up')
cli <- cli%>%
  select(bcr_followup_barcode,vital_status,
         days_to_death,days_to_last_followup)%>%
  distinct(bcr_followup_barcode,.keep_all = T)
# cli[is.na(cli)] <- 0

table(cli$vital_status)
dead_patient <- cli %>%
  dplyr::filter(vital_status == 'Dead')%>%
  dplyr::select(-days_to_last_followup) %>%
    rename(all_of(c(bcr_followup_barcode = 'Barcode',
             vital_status = 'fustat',
             days_to_death='futime'
              ))) %>%
  mutate(fustat=ifelse(fustat=='Dead',1,0))%>%
  mutate(futime=futime/365) 

#活的信息
alive_patient <-  cli %>%
  dplyr::filter(vital_status == 'Alive') %>%
  dplyr::select(-days_to_death) %>%
    rename(c(bcr_followup_barcode = 'Barcode',
             vital_status = 'fustat',
             days_to_last_followup='futime'
             )) %>%
  mutate(fustat=ifelse(fustat=='Dead',1,0))%>%
  mutate(futime=futime/365) 
#合并
survival_data <- rbind(dead_patient,alive_patient)
write.csv(survival_data,file="./survival/BLCA_surviv.csv")



BLCA_exp <- read.table("./survival/BLCA_exp.txt",sep = '\t',check.names = F,row.names = 1,header = T)
BLCA_exp <- data.frame(t(BLCA_exp))
rownames(BLCA_exp) <- gsub('\\.','\\-',rownames(BLCA_exp))

BLCA_surviv <- read.csv("./survival/BLCA_surviv.csv",row.names = 1) 

BLCA_exp$ID_exp <- substr(rownames(BLCA_exp), 1, 12)
BLCA_surviv$ID_sur <- substr(BLCA_surviv$Barcode, 1, 12)

map <- match(BLCA_exp$ID_exp,BLCA_surviv$ID_sur)

data3 <- BLCA_surviv[map,]##注意看这个
data4 <- cbind(BLCA_exp,data3)
data5 <- na.omit(data4)
desired_columns <- c("fustat","futime",'ID_exp','ID_sur')
all_columns <- c(desired_columns, setdiff(names(data5), desired_columns))
data6 <- data5[all_columns]

#############

# gene <- BLCA.table("../step1-TCSMP/myCytoscape/top10.txt",sep = "\t")
mydata <- data6[,colnames(data6) %in% c('fustat','futime','CHPF')]


######批量绘图
for (i in c('CHPF')) {
  ######根据均值将样本分为高低表达组
  mydata[,i] <- ifelse(mydata[,i]>=median(mydata[,i]),"high","low")
  splots <- list()
  km_fit <- survfit(Surv(futime,fustat)~mydata[,i],data = mydata)
  splots[[1]] <- ggsurvplot(km_fit,
                          xlab="Time_years",
                          ylab="Survival_probability",
                          pval = T,
                          conf.int = F,#置信区间设置
                          break.x.by=1,
                          legend.title=i,
                          legend.labs=c("high","low"),
                          palette = c("#A73030FF", "#0073C2FF"),
                          ggtheme = theme_bw(base_size = 18),#设置一种主题
                          risk.table.title="",#risk.table的标题为空
                          risk.table = TRUE, #显示risk.table
                          risk.table.col = "strata",#更改风险表的颜色
                          linetype = "strata", #更改线的类型
                          surv.median.line = "hv" 
                          )
  
  res <- arrange_ggsurvplots(splots,print = T,
                             ncol = 1,nrow = 1,risk.table.height = 0.25)
  ggsave(paste('./survival/',i,"All_surv.pdf",sep = "_"),res,width=7,height=7,dpi=600)#批量保存图片
}
class(splots)









