library(openxlsx)
library(venn) ### 7组及以下绘venn图
library(VennDiagram)###获得各个交集内的代谢物
library(UpSetR)###7组（不含）以上venn图变种
library(stringr)
library(readxl)






for (file in Sys.glob('../0 rawdata/*/*/TOTAL-Differentially Expressed Metabolites.xlsx')) {
  output<-strsplit(file,'/')[[1]]
  output<-paste0(output[3],'/',output[4])
  temp='.'
  for (ii in strsplit(output,'/')[[1]]) {
    temp<-paste0(temp,'/',ii)
    if(!dir.exists(temp)){dir.create(temp)}
  }
 
  datas=list()
  for (sheet in excel_sheets(file)[2:9]) {
    data0<-read.xlsx(file,sheet = sheet)
    data0<-data0[!is.na(data0$MS2.name),]
    # data0$SuperClass[is.na(data0$SuperClass)]<-'Others'
    data0<-data0[data0$MS2.score>0.6,]
    datas[[sheet]]<-data0
  }
  len<-max(unlist(lapply(datas, nrow)))+1
  data<-data.frame('a0'=rep('',len))
  for (ii in 1:length(datas)) {
    data[,paste0('a',ii)]<-c(datas[[ii]][,'MS2.name'],rep(NA,len-nrow(datas[[ii]])))
  }
  colnames(data)[-1]<-names(datas)
  data<-data[-len,-1]
  
  
  ###调用R包VennDiagram，获得各个交集内的代谢物
  sets<- get.venn.partitions(data)
  sets<-sets[sets['..count..']>0,]
  sets<-sets[!is.na(sets[ "..values.."]),]
  sets$..values..<-lapply(sets$..values..,function(i){i[!is.na(i)]})
  sets$..count..<-unlist(lapply(sets$..values..,length))
  sets$elements<-lapply(1:nrow(sets), function(i){
    paste(unlist(sets[i,'..values..']),collapse = ';;')
  }
  )
  sets<-sets[order(sets['..count..'],decreasing = TRUE),] 
  sets$'intersection'<-lapply(1:nrow(sets), function(i){
    ins<-strsplit(sets[i,'..set..'],')')[[1]][1]
    ins<-str_remove_all(ins,'\\(')#substring(ins,2,nchar(ins))
    gsub(pattern = '∩',replacement = '&',x=ins)
  }
  )
  insdata<-sets[,c('intersection','..count..',"elements")]
  colnames(insdata)<-c('intersection','count',"elements")
  insdata<-as.data.frame(lapply(insdata,as.character))
  write.csv(insdata,file = paste0(output,'/name_intersection.csv'),row.names = FALSE)
  write.xlsx(insdata,file = paste0(output,'/name_intersection.xlsx'),row.names = FALSE)
  insall<-sets[,!(colnames(sets) %in% c('..values..','..set..'))]
  insall<-as.data.frame(lapply(insall,as.character))
  write.csv(insall,file = paste0(output,'/intersection_all.csv'),row.names = FALSE)
  
  
  wuzhi<-unique(unlist(data))[!is.na(unique(unlist(data)))]
  groups<-colnames(data)
  venndata=matrix(nrow = length(wuzhi),ncol=ncol(data))
  rownames(venndata)<-wuzhi
  colnames(venndata)<-groups
  for (i in wuzhi){
    for(j in groups){
      if (i %in% data[,j]){
        venndata[i,j]<-1
      }else{
        venndata[i,j]<-0
      }
    }
  }
  venndata<-as.data.frame(venndata)
  write.csv(venndata,file = paste0(output,'/venn.csv'))
  
  # ### 
  # png(file = paste0(output,"/venn.png"),width = 1567*1.5,height = 1567,res = 300)
  # if(length(groups)<=4){
  #   venn(venndata,ilabels = TRUE, zcolor = "style",cexil = 1,cexsn=0.45,size=7,opacity=0.56)
  # }else if(length(groups)<=5){
  #   venn(venndata,ilabels = TRUE, zcolor = "style", ellipse = TRUE,opacity=0.56)
  # }else if(length(groups)<=7){
  #   venn(venndata,ilabels = TRUE, zcolor = "style", ellipse = FALSE,opacity=0.56)
  # }else{
  #   upset(venndata,sets.bar.color = "red",sets = groups,main.bar.color = "blue",order.by ='freq')
  # }
  # dev.off()
  # 
  # pdf(file = paste0(output,"/venn.pdf"))
  # if(length(groups)<=4){
  #   venn(venndata,ilabels = TRUE, zcolor = "style",opacity=0.56)
  # }else if(length(groups)<=5){
  #   venn(venndata,ilabels = TRUE, zcolor = "style", ellipse = TRUE,opacity=0.56)
  # }else if(length(groups)<=5){
  #   venn(venndata,ilabels = TRUE, zcolor = "style", ellipse = FALSE,opacity=0.56)
  # }else{
  #   upset(venndata,sets.bar.color = "red",sets = groups,main.bar.color = "blue",order.by ='freq')
  # }
  # dev.off()
  
}





