library(ggpubr)# 载入ggpubr这个包
#library(Hmisc) 

library(openxlsx)
library(xlsx)

setwd("E:\\project_Aurora\\内审\\谢会") #别忘记进入工作目录，就是你放文件的地方
# data <-read.xlsx(paste0("16s",k,".xlsx"), colNames = TRUE) 
# data_2 <- read.xlsx(paste0("metabolomics-FB-pos",k,".xlsx"), colNames = TRUE)
# aa <- t(data)
# bb <- t(data_2)
# aa_rownames = c("Group", sapply(rownames(aa), function(x){strsplit(x, "_")[[1]][1]}))
# aa_rownames <- aa_rownames[-2]
# aa_rownames <- as.data.frame(aa_rownames)
# aa_rownames
# 
# cc <- cbind(aa_rownames, aa, bb)
# cc1 <- cc[1,]
# cc2 <- as.matrix(cc1)
# colnames(cc) <- cc2
# cc <- cc[-1,]
# names(cc)
# names(cc) <- cc[1,]
# cc[1,]
# names(cc)
# cc <- cc[-1,]
#write.csv(cc, paste("scatter plot",k,".csv"))
data <-read.xlsx("01.xlsx",sheetIndex = 1,header = T)
#data1 <-read.xlsx("01.xlsx",sheetIndex = 1,header = F)#读取数据
#data <- as.matrix(data)
#lapply(data, as.numeric)
#group <- c("A",'A','A','B','B','B')
#data$group <- group
a <- names(data)
#b <- gsub("X","",a)
gene <- a[3:37]
meta <- a[38:41]


for(i in gene){
  for (j in meta){
    plotname <- paste(i, "-", j, sep = "")
    print(plotname)
    p <- ggscatter(data, x = i,y = j,shape = 20,size = 5,
                                 color = "Group",
                                 palette = c("#CD6090", "#1E90FF"),add.params = list(color = "black"),
                                 add = "reg.line",
                                 conf.int = T, cor.coef = T,show.legend.text = TRUE)+
      theme_bw() +
      theme(
        #修改字号
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 20 ,colour = "black"),
        #end
        legend.key = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 1)
      )
    # corr <- cor.test(data[,i],data[,j],method ="spearman")$estimate
    # Pvalue <- cor.test(data[,i],data[,j],method ="spearman")$p.value
    #  p <- p + annotate("text",x = max(data[,i])*0.8, y = max(data[,j])*1.05,label=paste0('rho = ', round(corr,3))) +
    #   annotate("text",x = max(data[,i])*0.8, y = max(data[,j])*1.02,label=paste0('p = ', round(Pvalue,3)))
  
    
  ggsave(paste(plotname,".jpg", sep = ""), p, width = 6, height = 7.5, units = 'in', dpi = 600)
  ggsave(paste(plotname,".pdf", sep = ""),p,width = 6,height = 7.5,units = 'in',dpi = 600)
  }
 }
 

  


