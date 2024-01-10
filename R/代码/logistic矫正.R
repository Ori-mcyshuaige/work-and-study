library(openxlsx)
library(stringr)
library(readxl)
library(reshape2)

meta<-read.xlsx("代谢.xlsx",rowNames = T,colNames = T,check.names = F,sep.names = " ")
pro<-read.xlsx("蛋白.xlsx",rowNames = T,colNames = T,check.names = F,sep.names = " ")
jzdata<-read.xlsx('临床矫正.xlsx',sheet = 2,rowNames = T)
jzdata$SEX<-factor(jzdata$SEX,levels = unique(jzdata$SEX))

ffdata<-data.frame(row.names = c("Estimate__(Intercept)","Estimate__pro","Estimate__sexM","Estimate__year",
                                 "Std. Error__(Intercept)","Std. Error__pro","Std. Error__sexM","Std. Error__year",
                                 "z value__(Intercept)","z value__pro","z value__sexM","z value__year",
                                 "Pr(>|z|)__(Intercept)","Pr(>|z|)__pro","Pr(>|z|)__sexM","Pr(>|z|)__year",
                                 "Qr(>|z|)__(Intercept)","Qr(>|z|)__pro","Qr(>|z|)__sexM","Qr(>|z|)__year"))

expdata<-pro[,rownames(jzdata)]
for(ii in rownames(expdata)){
  data<-data.frame(pro=as.numeric(expdata[ii,]),
                   sex=jzdata[,'SEX'],
                   year=jzdata[,'year'],
                   # BMI=jzdata[,"BMI"],
                   group=jzdata[,'G'],
                   row.names = rownames(jzdata))
  data$group<-factor(data$group,levels = unique(data$group))
  lghg<-glm(group~pro+sex+year,family = binomial(link = "logit"),data = data,maxit=200)
  idata<-data.frame(coef(summary(lghg)),check.names = F)
  idata$`Qr(>|z|)`<-p.adjust(idata$`Pr(>|z|)`,method = 'BH')
  iidata<-melt(as.matrix(idata))
  rownames(iidata)<-paste0(iidata$Var2,'__',iidata$Var1)
  # iidata<-iidata[rownames(ffdata),]
  ffdata[,ii]<-iidata[rownames(ffdata),'value']
}
fffdata<-data.frame(t(ffdata),check.names = F)
write.xlsx(fffdata,'蛋白年龄性别_adjust.xlsx',rowNames=T)

# diffdata0<-diffexp[rownames(fffdata)[fffdata$`Qr(>|z|)__pro`<0.05],]
# write.xlsx(diffdata0,paste0(sheet,'.xlsx'),rowNames=T)

diffdata<-expdata[rownames(fffdata)[fffdata$`Qr(>|z|)__pro`<0.05],]
write.xlsx(diffdata,'蛋白年龄性别_exp.xlsx',rowNames=T)
