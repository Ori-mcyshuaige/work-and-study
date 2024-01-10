library(openxlsx)
# library(XLConnect)
setwd("D:/马驰宇/181/戴主任/工作/生信流程/项目/分析/AQ-MRZ20200817-FC")
mean<-read.xlsx("POS-Mean.xlsx",sheet = 2,colNames = T,rowNames = T,check.names = F,sep.names = " ")
mean<-mean[!is.na(mean[,1]),]
rownames(mean)<-mean[,1]
group<-c("Y_2019 vs Y_2018","Y_2016 vs Y_2015","Y_2018 vs Y_2017","Y_2016 vs D_2016","Y_2015 vs Y_2014","Y_2017 vs Y_2016","Y_2016 vs C_2016","D_2016 vs C_2016")
wb <-createWorkbook()
modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
hs <- createStyle(
  fontSize = 10, fontColour = 'white', fgFill = '#008CCE', 
  valign = 'center', halign = 'center'
)
for(i in 2:9){
  deg<-read.xlsx("POS-Differentially Expressed Metabolites.xlsx",sheet = i,colNames = T,rowNames = T,check.names = F,sep.names = " ")
  deg<-deg[!is.na(deg[,1]),1]
  gname<-strsplit(group[i-1]," vs ",fixed = T)
  gn1<-gname[[1]][1]
  gn2<-gname[[1]][2]
  deg.mean<-mean[deg,c(1,grep(gn1,colnames(mean)),grep(gn1,colnames(mean)))]
  addWorksheet(
    wb, sheetName = group[i-1],
    gridLines = F, tabColour = '#008CCE'
  )
  writeDataTable(
    wb, sheet = i-1, deg.mean, headerStyle = hs,
    tableStyle = 'TableStyleMedium27', withFilter = F
  )
  setColWidths(wb, sheet = i-1, cols = 1:ncol(deg.mean), 'auto')
  setRowHeights(wb, sheet = i-1, rows = 1, 24)
  setRowHeights(wb, sheet = i-1, rows = 2:(nrow(deg.mean) + 1), 16)
  addStyle(
    wb, sheet = i-1, rows = 2:(nrow(deg.mean) + 1),
    cols = 1:ncol(deg.mean), gridExpand = T, stack = T,
    style = createStyle(valign = 'center')
  )
}
saveWorkbook(wb,overwrite = T,"hierarchical clustering data matrix.xlsx")
  
  