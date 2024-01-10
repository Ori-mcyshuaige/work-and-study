library(openxlsx)
library(venn) ### 7组及以下绘venn图
library(VennDiagram)###获得各个交集内的代谢物
library(UpSetR)###7组（不含）以上venn图变种
library(stringr)
library(readxl)
library(grDevices)

# excle<-list.files()
# data<-data.frame(matrix(nrow = 1300,ncol = 7))
# colnames(data)<-excle
# for(i in excle){
#   raw<-read.xlsx(i,sheet = 5,rowNames = T,sep.names = " ")
#   select<-NULL
#   for(j in colnames(raw)){
#     p<-rownames(raw)[raw[,j]<0.05]
#     select<-c(select,p)
#   }
#   select<-unique(select)
#   data[,i][1:length(select)]<-select
# }
# colnames(data)<-unlist(strsplit(excle,".xlsx"))
# setwd("D:/马驰宇/181/戴主任/工作/生信流程/项目/分析/何院专辑/某个项目1")
data<-read.xlsx("01.节律与差异节律/venndata.xlsx",sheet = 1,startRow = 1,colNames = TRUE,rowNames = FALSE,check.names = FALSE,sep.names = " ")
venn<-vector(mode = "list")
for(i in 1:length(colnames(data))){
  venn[[colnames(data)[i]]]=data[,i][!is.na(data[,i])]
}
sets<- get.venn.partitions(venn)
sets<-sets[sets['..count..']>0,]
sets$elements<-lapply(1:nrow(sets), function(i){
  paste(unlist(sets[i,'..values..']),collapse = ';')
}
)
sets$'intersection'<-lapply(1:nrow(sets), function(i){
  ins<-strsplit(sets[i,'..set..'],')')[[1]][1]
  ins<-str_remove_all(ins,'\\(')#substring(ins,2,nchar(ins))
  gsub(pattern = '∩',replacement = '&',x=ins)
}
)
insdata<-sets[,c('intersection','..count..',"elements")]
colnames(insdata)<-c('intersection','count',"elements")
insdata<-as.data.frame(lapply(insdata,as.character))
png(file = "venn2.png",width = 3100*1.5,height = 3100,res = 600)
venn(venn,zcolor="style",opacity = 0.3, plotsize = 15, ilcs = 1, sncs = 1.2,
     borders = T, box = FALSE, par = T)
dev.off()
pdf(file = "venn2.pdf",width = 6.2*1.5,height = 6.2)
venn(venn,zcolor="style",opacity = 0.3, plotsize = 15, ilcs = 1, sncs = 1.2,
     borders = T, box = FALSE, par = T)
dev.off()
write.csv(insdata,"venn result.csv")

venn1<-venn.diagram(venn,"venn1.png",fill=jet(2),alpha=0.5,cat.fontfamily = "sans",fontfamily = "sans",margin = 0.05,
                    col = "transparent",lwd = 3,cex = 1.2)
                    
pdf(file="venn1.pdf")
venn1<-venn.diagram(venn,filename = NULL,fill=jet(2),alpha=0.5,cat.fontfamily = "sans",fontfamily = "sans",margin = 0.05,
                    col = "transparent",lwd = 3,cex = 1.2)
grid.draw(venn1)
dev.off()



# # --多组venn图
# data<-read.xlsx("Venn1.xlsx",sheet = 2,startRow = 1,colNames = TRUE,rowNames = T,check.names = FALSE,sep.names = " ")
# venn<-vector(mode = "list")
# for(i in 1:length(colnames(data))){
#   venn[[colnames(data)[i]]]=data[,i][!is.na(data[,i])]
# }
# sets<- get.venn.partitions(venn)
# sets<-sets[sets['..count..']>0,]
# sets$elements<-lapply(1:nrow(sets), function(i){
#   paste(unlist(sets[i,'..values..']),collapse = ';;')
# }
# )
# sets$'intersection'<-lapply(1:nrow(sets), function(i){
#   ins<-strsplit(sets[i,'..set..'],')')[[1]][1]
#   ins<-str_remove_all(ins,'\\(')#substring(ins,2,nchar(ins))
#   gsub(pattern = '∩',replacement = '&',x=ins)
# }
# )
# insdata<-sets[,c('intersection','..count..',"elements")]
# colnames(insdata)<-c('intersection','count',"elements")
# insdata<-as.data.frame(lapply(insdata,as.character))
# write.csv(insdata,"venn result.csv")
# png(file = "venn.png",width = 3100*1.5,height = 3100,res = 600,family = "A")
# upset(data,nsets = ncol(data),nintersects = NA,set_size.show = T,mainbar.y.label = NULL,sets.x.label = NULL,set_size.numbers_size = 8)
# dev.off()

