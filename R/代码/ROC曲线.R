library(pROC)
library(dplyr)
library(pals)
library(openxlsx)

###标点  print.thres="best",thresholds="best",

data<-read.xlsx("ROC.xlsx",sheet = 1,check.names = F,sep.names = " ",rowNames = F)

##组合ROC
png("ROC.png",width = 10,height = 10,units = "in",res = 600)
# par(xpd=T,pin=c(4,4),mar=c(6,5,5,8)+0.1)
p1<-roc(data$group,data[,high[1]],percent = T,ci=T,plot=T,grid=F,print.auc=F,col=color[1],legacy.axes=T)
p2<-roc(data$group,data[,high[2]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[2],legacy.axes=T)
p3<-roc(data$group,data[,high[3]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[3],legacy.axes=T)
p4<-roc(data$group,data[,high[4]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[4],legacy.axes=T)
p5<-roc(data$group,data[,high[5]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[5],legacy.axes=T)
p6<-roc(data$group,data[,high[6]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[6],legacy.axes=T)
p7<-roc(data$group,data[,high[7]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[7],legacy.axes=T)
p8<-roc(data$group,data[,high[8]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[8],legacy.axes=T)
p9<-roc(data$group,data[,high[9]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[9],legacy.axes=T)
p10<-roc(data$group,data[,high[10]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[10],legacy.axes=T)
p11<-roc(data$group,data[,high[11]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[11],legacy.axes=T)
p12<-roc(data$group,data[,high[12]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[12],legacy.axes=T)
p13<-roc(data$group,data[,high[13]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[13],legacy.axes=T)
p14<-roc(data$group,data[,high[14]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[14],legacy.axes=T)
p15<-roc(data$group,data[,high[15]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[15],legacy.axes=T)
p16<-roc(data$group,data[,high[16]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[16],legacy.axes=T)
p17<-roc(data$group,data[,high[17]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[17],legacy.axes=T)
# xy=par("usr")
legend('bottomright',
       # x=par("usr")[2],y=par("usr")[4],
       title = '',
       legend = c(paste0(colnames(data)[high[1]],":",sprintf("%.2f (%.2f-%.2f) %%", p1$auc, p1$ci[1],p1$ci[3])),
                  paste0(colnames(data)[high[2]],":",sprintf("%.2f (%.2f-%.2f) %%", p2$auc, p2$ci[1],p2$ci[3])),
                  paste0(colnames(data)[high[3]],":",sprintf("%.2f (%.2f-%.2f) %%", p3$auc, p3$ci[1],p3$ci[3])),
                  paste0(colnames(data)[high[4]],":",sprintf("%.2f (%.2f-%.2f) %%", p4$auc, p4$ci[1],p4$ci[3])),
                  paste0(colnames(data)[high[5]],":",sprintf("%.2f (%.2f-%.2f) %%", p5$auc, p5$ci[1],p5$ci[3])),
                  paste0(colnames(data)[high[6]],":",sprintf("%.2f (%.2f-%.2f) %%", p6$auc, p6$ci[1],p6$ci[3])),
                  paste0(colnames(data)[high[7]],":",sprintf("%.2f (%.2f-%.2f) %%", p7$auc, p7$ci[1],p7$ci[3])),
                  paste0(colnames(data)[high[8]],":",sprintf("%.2f (%.2f-%.2f) %%", p8$auc, p8$ci[1],p8$ci[3])),
                  paste0(colnames(data)[high[9]],":",sprintf("%.2f (%.2f-%.2f) %%", p9$auc, p9$ci[1],p9$ci[3])),
                  paste0(colnames(data)[high[10]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[11]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[12]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[13]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[14]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[15]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[16]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[17]],":",sprintf("%.2f (%.2f-%.2f) %%", p11$auc, p11$ci[1],p11$ci[3]))
       ),
       lty = 1,
       lwd = 2,
       col=color,#'tomato','orange',"purple","#708090","brown","#3CB371","#54FF9F","#6495ED","#191970","#FFC1C1","#8B658B"
       bty = 'n',##图例框：o,n
)
dev.off()

pdf("ROC.pdf",width = 10,height = 10)
# par(xpd=T,pin=c(4,4),mar=c(6,5,5,8)+0.1)
p1<-roc(data$group,data[,high[1]],percent = T,ci=T,plot=T,grid=F,print.auc=F,col=color[1],legacy.axes=T)
p2<-roc(data$group,data[,high[2]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[2],legacy.axes=T)
p3<-roc(data$group,data[,high[3]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[3],legacy.axes=T)
p4<-roc(data$group,data[,high[4]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[4],legacy.axes=T)
p5<-roc(data$group,data[,high[5]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[5],legacy.axes=T)
p6<-roc(data$group,data[,high[6]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[6],legacy.axes=T)
p7<-roc(data$group,data[,high[7]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[7],legacy.axes=T)
p8<-roc(data$group,data[,high[8]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[8],legacy.axes=T)
p9<-roc(data$group,data[,high[9]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[9],legacy.axes=T)
p10<-roc(data$group,data[,high[10]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[10],legacy.axes=T)
p11<-roc(data$group,data[,high[11]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[11],legacy.axes=T)
p12<-roc(data$group,data[,high[12]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[12],legacy.axes=T)
p13<-roc(data$group,data[,high[13]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[13],legacy.axes=T)
p14<-roc(data$group,data[,high[14]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[14],legacy.axes=T)
p15<-roc(data$group,data[,high[15]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[15],legacy.axes=T)
p16<-roc(data$group,data[,high[16]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[16],legacy.axes=T)
p17<-roc(data$group,data[,high[17]],percent = T,ci=T,plot=T,grid=F,print.auc=F,add=T,col=color[17],legacy.axes=T)
# xy=par("usr")
legend('bottomright',
       # x=par("usr")[2],y=par("usr")[4],
       title = '',
       legend = c(paste0(colnames(data)[high[1]],":",sprintf("%.2f (%.2f-%.2f) %%", p1$auc, p1$ci[1],p1$ci[3])),
                  paste0(colnames(data)[high[2]],":",sprintf("%.2f (%.2f-%.2f) %%", p2$auc, p2$ci[1],p2$ci[3])),
                  paste0(colnames(data)[high[3]],":",sprintf("%.2f (%.2f-%.2f) %%", p3$auc, p3$ci[1],p3$ci[3])),
                  paste0(colnames(data)[high[4]],":",sprintf("%.2f (%.2f-%.2f) %%", p4$auc, p4$ci[1],p4$ci[3])),
                  paste0(colnames(data)[high[5]],":",sprintf("%.2f (%.2f-%.2f) %%", p5$auc, p5$ci[1],p5$ci[3])),
                  paste0(colnames(data)[high[6]],":",sprintf("%.2f (%.2f-%.2f) %%", p6$auc, p6$ci[1],p6$ci[3])),
                  paste0(colnames(data)[high[7]],":",sprintf("%.2f (%.2f-%.2f) %%", p7$auc, p7$ci[1],p7$ci[3])),
                  paste0(colnames(data)[high[8]],":",sprintf("%.2f (%.2f-%.2f) %%", p8$auc, p8$ci[1],p8$ci[3])),
                  paste0(colnames(data)[high[9]],":",sprintf("%.2f (%.2f-%.2f) %%", p9$auc, p9$ci[1],p9$ci[3])),
                  paste0(colnames(data)[high[10]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[11]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[12]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[13]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[14]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[15]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[16]],":",sprintf("%.2f (%.2f-%.2f) %%", p10$auc, p10$ci[1],p10$ci[3])),
                  paste0(colnames(data)[high[17]],":",sprintf("%.2f (%.2f-%.2f) %%", p11$auc, p11$ci[1],p11$ci[3]))
       ),
       lty = 1,
       lwd = 2,
       col=color,#'tomato','orange',"purple","#708090","brown","#3CB371","#54FF9F","#6495ED","#191970","#FFC1C1","#8B658B"
       bty = 'n',##图例框：o,n
)
dev.off()


###单独ROC
for(i in 3:length(colnames(data))){
  png(paste0(colnames(data)[i],".png"),width = 7,height = 7,units = "in",res = 600)
  p<-roc(data$group,data[,i],ci=T,main = colnames(data)[i],plot=T,grid=F,print.auc=T,col="tomato",legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
  pdf(paste0(colnames(data)[i],".pdf"),width = 7,height = 7)
  p<-roc(data$group,data[,i],ci=T,main = colnames(data)[i],plot=T,grid=F,print.auc=T,col="tomato",legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
}

filelist <- c('risk_GSE32280.csv','risk_GSE76826.csv','risk_GSE98793.csv')
for(i in filelist){
  data <- read.csv(i,row.names = 1)
  png(paste0(strsplit(i,'\\.')[[1]][1],".png"),width = 7,height = 7,units = "in",res = 600)
  p<-roc(data$group,data$risk_score,ci=T,main = colnames(data)[i],plot=T,grid=F,print.auc=T,col="red",legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
  pdf(paste0(strsplit(i,'\\.')[[1]][1],".pdf"),width = 7,height = 7)
  p<-roc(data$group,data$risk_score,ci=T,main = colnames(data)[i],plot=T,grid=F,print.auc=T,col="red",legacy.axes=T,print.thres="best",thresholds="best")
  dev.off()
}
