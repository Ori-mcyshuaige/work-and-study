###联合ROC
# setwd("C:/Users/小马很酷/Desktop/新建文件夹 (2)/周学敏")
library(openxlsx)
library(pROC)

data<-read.xlsx("IgAN三个蛋白表达量及分组.xlsx",sheet = 2,colNames = T,check.names = F)###原始数据
data$Group<-as.factor(data$Group)
met<-names(data)[-(1:2)]
met1<-sapply(met, function(x){x<-paste0("data$","`",x,"`")})
names<-gsub("\\."," ",met)
combine<-function(m){     #获取随机组合列表
  c<-list()
  for (i in 3:length(m)){
    b<<-combn(m,i)
    c[[i]]<-apply(b, 2 ,function(x){x<<-paste0(x,collapse = "+")})
  }
  return(unlist(c))
}
# 1:length(m)
com.names<-combine(names)
com.names1<- gsub(
  "/|:|\\*|\\?|\"|<|>|\\|", "", com.names, perl = T
)
com.met<-combine(met1)
matrix.auc<-matrix(0,ncol=2,nrow=length(com.met))
for (j in 1:length(com.names1)){
  f <- as.formula(paste0("data$Group~", com.met[j]))
  g<-glm(f,family=binomial(link='logit'))###逻辑回归
  pre<-predict(g)
  pre<-exp(pre)/(1+exp(pre))
  p2<-roc(data$Group~pre,
         smoothed = TRUE,
         ci=TRUE, 
         plot=T, 
         max.auc.polygon=TRUE, 
         grid=F,
         legacy.axes=T,print.thres="best",thresholds="best",
         print.auc=TRUE
  )
  Roc = paste0("ROC","_",com.names1[j],".png")
  Roc1 = paste0("ROC","_",com.names1[j],".pdf")
  png(Roc,res = 300,width = 2000,height = 2000)
  roc_plot<-plot(p2,
                 ci=TRUE,
                 max.auc.polygon=F, 
                 grid=F,
                 print.auc=TRUE,print.thres="best",thresholds="best",
                 legacy.axes=T,
                 lwd=2
  )
  dev.off()

  pdf(Roc1,width = 5,height = 5)
  roc_plot<-plot(p2,
                 ci=TRUE,
                 max.auc.polygon=F, 
                 grid=F,
                 print.auc=TRUE,print.thres="best",thresholds="best",
                 legacy.axes=T,
                 lwd=2
  )
  dev.off()
  matrix.auc[j,1]<-round(as.numeric(p2$ci)[2],3)
  matrix.auc[j,2]<-paste0(round(as.numeric(p2$ci)[1],3),"--",round(as.numeric(p2$ci)[3],3))
}
# colnames(matrix.auc)<-c("AUC","95%CI")
# row.names(matrix.auc)<-com.names
# matrix.auc<-as.data.frame(matrix.auc)
write.xlsx(matrix.auc,"AUC.xlsx",rowNames=T,colNames=T)


# png("total图3.png",res = 300,width = 2000,height = 2000)
# 
# roc_plot<-plot(p1,
#                ci=TRUE,
#                max.auc.polygon=TRUE, 
#                grid=TRUE,
#                col="red",
#                lwd=2,
#                bg="red"
# )
# lines(p2,col="yellow")
# lines(p3,col="blue")
# lines(p4,col="pink")
# lines(p5,col="green")
# lines(p6,col="orange")
# lines(p7,col="purple")
# legend(0.55,0.23,c("Pyruvic acid  AUC:0.888(0.800-0.976)",
#                  "succinate semialdehyde 2  AUC:0.853(0.735-0.972)",
#                  "fumaric acid  AUC:0.717(0.561-0.873)",
#                  "L-Malic acid  AUC:0.808(0.690-0.937)",
#                  "alpha-ketoglutaric acid  AUC:0.812(0.687-0.937)",
#                  "citric acid  AUC:0.948(0.898-0.998)",
#                  "Total  AUC:0.738(0.607-0.869)"
#                  ),col=c("red","yellow","blue","pink","green","orange","purple"),cex=0.6,lty = 1)
# text(0.6,0.6,"AUC=0.888",col="red",cex = 0.75)
# text(0.6,0.5,"AUC=0.853",col="yellow",cex = 0.75)
# text(0.6,0.4,"AUC=0.717",col="blue",cex = 0.75)
# text(0.6,0.3,"AUC=0.808",col="pink",cex = 0.75)
# text(0.2,0.8,"AUC=0.812",col="green",cex = 0.75)
# text(0.2,0.6,"AUC=0.948",col="orange",cex = 0.75)
# text(0.2,0.4,"AUC=0.738",col="purple",cex = 0.75)
# dev.off()
