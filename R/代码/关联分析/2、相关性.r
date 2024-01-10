#相关性计算脚本(会出3个结果，P值，Cor值，和omics)
library(openxlsx)
# setwd('C:\\Users\\Administrator\\Desktop\\关联分析\\相关性分析')
options(stringsAsFactors = F)
#workDir1 <- "C:\\Users\\Administrator\\Desktop\\关联分析\\相关性分析\\"
#workDir2 <- "C:\\Users\\Administrator\\Desktop\\关联分析\\相关性分析\\"
sampleSet1 <- read.xlsx("Metabonomics.xlsx",colNames = T,rowNames = T,sheet = 3)
sampleSet2 <- read.xlsx("16s.xlsx",colNames = T,rowNames = T,sheet = 3)

n1 <- nrow(sampleSet1)
n2 <- nrow(sampleSet2)
node1 <- gsub("\\."," ",rownames(sampleSet1))
node2 <- gsub("\\."," ",rownames(sampleSet2))
corrList <- list()
m <- matrix(nrow = n1,ncol = n2)
dfM <- data.frame(m)
rownames(dfM) <- rownames(sampleSet1)
colnames(dfM) <- rownames(sampleSet2)
matrixCorr <- dfM
matrixPvalue <- dfM
for(i in 1:n1){
  a<-i
  for(j in 1:n2){
    x <- as.numeric(unlist(sampleSet1[i,]))
    y <- as.numeric(unlist(sampleSet2[j,]))
    corr <- cor.test(x,y,method = "spearman")
    matrixCorr[i,j] <- corr$estimate
    matrixPvalue[i,j] <- corr$p.value
    element <- c(node1[i],node2[j],as.numeric(corr$estimate),as.numeric(corr$p.value))
    corrList <- c(corrList,list(element))
  }
}

net <- t(data.frame(corrList))
colnames(net) <- c("node1","node2","corr","corr.p")
write.xlsx(net,"network.xlsx",col.names = T,row.names = F)
write.xlsx(matrixCorr,"matrixCorr.xlsx",col.names = T,row.names = T)
write.xlsx(matrixPvalue,"matrixPvalue.xlsx",col.names = T,row.names = T)

nodeAttr <- data.frame(c(node1,node2),"omics")
colnames(nodeAttr) <- c("node","omics")
nodeAttr[nodeAttr$node %in% node1,"omics"] <- "Metabonomics"
nodeAttr[nodeAttr$node %in% node2,"omics"] <- "Transcriptomics"
write.xlsx(nodeAttr,"nodeAttr.xlsx",col.names = T,row.names = F)
##############绘图
library(ggplot2)

    corr <- read.xlsx("network.xlsx",colNames = T,rowNames = F)
    corr$corr <- as.numeric(corr$corr)
    corr$corr.p <- as.numeric(corr$corr.p)
    p <- ggplot(corr, aes(node2, node1)) + geom_tile(aes(fill = corr), color = "white", size = 0.1) + scale_fill_gradient2(
      breaks = c(-1, -0.5, 0, 0.5, 1),
      limits = c(-1, 1),
      labels = c(-1, -0.5, 0, 0.5, 1),
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar"
    ) +
      annotate(
        "text",
        x = corr$node2[corr$corr.p < 0.05] ,
        y = corr$node1[corr$corr.p < 0.05],
        label = "*",
        size = 3
      ) +
      theme(
        axis.text.x = element_text(
          angle = -90,
          hjust  = 0,
          vjust = 0.5
        ),
        axis.text = element_text(colour = "#000000"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA),
        legend.justification =  c("right","top")
      ) +
      labs(x = NULL ,y = NULL)
    picWidth <- length(unique(corr$node2))/4 + max(nchar(corr$node1))/20+2
    picHeight <- length(unique(corr$node1))/4 + max(nchar(corr$node2))/20
    ggsave(
      "corr.pdf",
      plot = p,
      width = picWidth,
      height = picHeight,
      limitsize = FALSE
    )
   ggsave(
     "corr.tiff",
     plot = p,
     width = picWidth,
     height = picHeight,
     limitsize = FALSE
   )####图片信息太多会报错

  