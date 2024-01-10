library(openxlsx)
library(readxl)
library(psych)
library(FSA)
library(vegan)
library(stringr)
options(stringsAsFactors = F,na.rm=T)


before<-paste0('Before',seq(1,63,by = 2))
after<-paste0('After',seq(2,64,by = 2))



datatest0<-t(data.frame(Statistics=c('Statistics','Statistics','Mean','Mean','wilcox Test','wilcox Test','PERMANOVA','PERMANOVA'),
                         row.names = c('Class','Characteristics','Before','After','Pvalue(After-Vs-Before)','After-Before','R2','P-value')))
datatest0<-as.data.frame(datatest0)


data1<-read.xlsx('临床指标汇总.xlsx',rowNames=T)
datatest1<-datatest0[rep(1,ncol(data1)),]
datatest1[,1]<-'临床指标'
datatest1[,2]<-colnames(data1)
datatest1[,3]<-paste0(
  mapply(mean, matrix(data1[before,]),na.rm=T),
  '(',
  mapply(min, matrix(data1[before,]),na.rm=T),
  '~',
  mapply(max, matrix(data1[before,]),na.rm=T),
  ')'
)
datatest1[,4]<-paste0(
  mapply(mean, matrix(data1[after,]),na.rm=T),
  '(',
  mapply(min, matrix(data1[after,]),na.rm=T),
  '~',
  mapply(max, matrix(data1[after,]),na.rm=T),
  ')'
)
datatest1[,6]<-paste0(
  mapply(mean, matrix(data1[after,]-data1[before,]),na.rm=T),
  '(',
  mapply(min, matrix(data1[after,]-data1[before,]),na.rm=T),
  '~',
  mapply(max, matrix(data1[after,]-data1[before,]),na.rm=T),
  ')'
)

datatest1[,5]<-unlist(sapply(colnames(data1), function(i){
  wilcox.test(data1[before,i],data1[after,i],paired = T)$p.value
}))



data2<-data.frame(Group=unlist(sapply(rownames(data1), function(i){paste0(str_extract_all(i,'\\D')[[1]],collapse = '')})),
                  row.names = rownames(data1))

for (ii in colnames(data1)) {
  if(length(which(is.na(data1[,ii])))==0){
    ad<-adonis(data1[,ii]~Group,data2, permutations = 999,na.rm=T)
  }else{
    adonis(data1[which(!is.na(data1[,ii])),ii]~Group,
           data.frame(Group=data2[which(!is.na(data1[,ii])),],row.names=rownames(data2)[which(!is.na(data1[,ii]))]), 
           permutations = 999,na.rm=T)
  }
  
  datatest1[datatest1$Characteristics==ii,7]<-ad$aov.tab[1,5]
  datatest1[datatest1$Characteristics==ii,8]<-ad$aov.tab[1,6]
}
datatest<-rbind(datatest0,datatest1)




# datatest<-datatest0
for (sheet in excel_sheets('体检信息汇总.xlsx')) {
  data1<-read.xlsx('体检信息汇总.xlsx',sheet=sheet,rowNames=T)
  data1[data1=='']<-NA
  ###秩和测验
  ###perno
  datatest1<-datatest0[rep(1,ncol(data1)),]
  datatest1[,1]<-sheet
  datatest1[,2]<-colnames(data1)
  datatest1[,3]<-paste0(
    mapply(mean, matrix(data1[before,]),na.rm=T),
    '(',
    mapply(min, matrix(data1[before,]),na.rm=T),
    '~',
    mapply(max, matrix(data1[before,]),na.rm=T),
    ')'
  )
  datatest1[,4]<-paste0(
    mapply(mean, matrix(data1[after,]),na.rm=T),
    '(',
    mapply(min, matrix(data1[after,]),na.rm=T),
    '~',
    mapply(max, matrix(data1[after,]),na.rm=T),
    ')'
  )
  datatest1[,6]<-paste0(
    mapply(mean, matrix(data1[after,]-data1[before,]),na.rm=T),
    '(',
    mapply(min, matrix(data1[after,]-data1[before,]),na.rm=T),
    '~',
    mapply(max, matrix(data1[after,]-data1[before,]),na.rm=T),
    ')'
  )
  
  datatest1[,5]<-unlist(sapply(colnames(data1), function(i){
    wilcox.test(data1[before,i],data1[after,i],paired = T)$p.value
  }))
  
  
  
  data2<-data.frame(Group=unlist(sapply(rownames(data1), function(i){paste0(str_extract_all(i,'\\D')[[1]],collapse = '')})),
                    row.names = rownames(data1))
  
  for (ii in colnames(data1)) {
    if(length(which(is.na(data1[,ii])))==0){
      ad<-adonis(data1[,ii]~Group,data2, permutations = 999,na.rm=T)
    }else{
      adonis(data1[which(!is.na(data1[,ii])),ii]~Group,
             data.frame(Group=data2[which(!is.na(data1[,ii])),],row.names=rownames(data2)[which(!is.na(data1[,ii]))]), 
             permutations = 999,na.rm=T)
    }
    
    datatest1[datatest1$Characteristics==ii,7]<-ad$aov.tab[1,5]
    datatest1[datatest1$Characteristics==ii,8]<-ad$aov.tab[1,6]
  }
  datatest<-rbind(datatest,datatest1)
  
}





## Create a new workbook 
wb <- createWorkbook()
## Add a worksheet 
hs <- createStyle(
  fontSize = 20, fontColour = 'white', fgFill = '#0064FF',borderStyle='double',
  valign = 'center', halign = 'center',border="TopBottomLeftRight", borderColour = "black"
)
addWorksheet(wb, "Statistics", gridLines = F, tabColour = '#0064FF') 

writeDataTable(
  wb, sheet = "Statistics", datatest, headerStyle = hs,startRow = 1,
  tableStyle = "TableStyleMedium25", withFilter = F
)



for (jj in unique(datatest[1,])[-1]) {
  mergeCells(wb, "Statistics", cols =which(datatest[1,]==jj), rows = 2)  
}

for (mm in unique(datatest$Class[-1])) {
  mergeCells(wb, "Statistics", cols = 1, rows = 1+which(datatest$Class==mm))  
}

headerStyle1 <- createStyle(fontSize = 14, fontColour = "#FFFFFF",  valign = 'center', halign = 'center', fgFill = "#4F81BD",
                            border="TopBottomLeftRight", borderColour = "black")
addStyle(wb, sheet = 'Statistics', headerStyle1, rows = 1:nrow(datatest)+1, cols = 1, gridExpand = TRUE)
headerStyle2 <- createStyle(fontSize = 9, fontColour = "black", halign = "center",  valign = 'center', 
                            border="TopBottomLeftRight", borderColour = "black",borderStyle='dashed')
addStyle(wb, sheet = 'Statistics', headerStyle2, rows = 1:nrow(datatest)+2, cols = 2:ncol(datatest), gridExpand = TRUE)
setColWidths(wb, sheet = "Statistics", cols = 1:ncol(datatest), 'auto')
setRowHeights(wb, sheet = "Statistics", rows = 1,24)
setRowHeights(wb, sheet = "Statistics", rows = 1:nrow(datatest)+1, 16)
saveWorkbook(wb,paste0('data.xlsx'),overwrite = TRUE)