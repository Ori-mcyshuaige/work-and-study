indir<-"D:/马驰宇/181/戴主任/工作/生信流程/项目/分析/杨勇 AQ-YY20200929-QTOF-XSXZ/POS/CHOW-GFP-CHOW-OE/kegg"
hitMin <- 1
library(rvest)
library(openxlsx)
library(tidyverse)
#不同浏览器html pathway分隔符 (目前支持360，火狐，谷歌）
{
  sepP1<-c("\\)\n\n","\\)\r\n\r\n")
  sepCG1<-c("[a-z]\n","[a-z]\r\n")
}
cat("start to deal with the kegg hit list,keggHitMin is set to",hitMin,"\n")
file<-list.files(path = indir,pattern = "htm")
fileTest<-grep("\\.html",file)
{
  if (!is_empty(fileTest)){
    oriCode<-read_html(paste0(indir,"/KEGG Mapper Search Result.html"))
}
  else {
    oriCode<-read_html(paste0(indir,"/KEGG Mapper Search Result.htm"))
  }
}
nodeslink<-html_nodes(oriCode,"li")
No<-length(nodeslink)
sepTest<-function(data,sepP0=sepP1,sepCG0=sepCG1){
  sepPath<-sapply(sepP0,function(i){
    pk<-grepl(i,data)
    if (is_empty(pk)){
      pk<-NA
    }
    else {pk<-pk}
    i<-pk
  })
  sepCpdGene<-sapply(sepCG0,function(i){
    pk<-grepl(i,data)
    if (is_empty(pk)){
      pk<-NA
    }
    else {pk<-pk}
    i<-pk
  })
  if (any(sepPath)){
    sepP<-sepP0[which(sepPath==TRUE)]
    sepP<-sub("\\\\)","",sepP)
  }
  else {sepP<-NA}
  if (any(sepCpdGene)){
    sepCG<-sepCG0[which(sepCpdGene==TRUE)]
    sepCG<-sub("\\[a-z\\]","",sepCG)
  }
  else {sepCG<-NA}
  return(list(sepP=sepP,sepCG=sepCG))
}
sepfilter<-sepTest(data = html_text(nodeslink[5]))
if (any(is.na(sepfilter))){
  stop("wrong separation character of pathway information, please check separation character")
}
pathwayInfo<-lapply(nodeslink[4:No],function(i){
  pathName<-unlist(strsplit(html_text(i),sepfilter$sepP))[1]
  cpdgeneInfo<-unlist(strsplit(html_text(i),sepfilter$sepP))[2]
  cpdgene<-sapply(cpdgeneInfo,function(j){
    j<-unlist(strsplit(j,sepfilter$sepCG))
  })
  cpdgene<-gsub("^ *","",cpdgene)
  pathHits<-paste0(cpdgene,collapse = "||")
  oriBaits<-gsub(" .*$","",cpdgene)
  oriBaits<-gsub("cpd:","",oriBaits)
  oriBaits<-gsub("ko:","",oriBaits)
  oriBaits<-paste0(oriBaits,collapse = "|")
  return(list(pathName=pathName,pathHits=pathHits,oriBaits=oriBaits))
})
pathway<-data.frame(
  "pathwayName"=sapply(pathwayInfo,`[[`,1,USE.NAMES = F),
  "pathwayHits"=sapply(pathwayInfo,`[[`,2,USE.NAMES = F),
  "oriBaits"=sapply(pathwayInfo,`[[`,3,USE.NAMES = F)
)
pathway_No<-sapply(pathway$pathwayName,function(i){
  name<-sub("^.*\\(","",i)
  count<-sub("\\)","",name)
  i<-count
})
pathway<-pathway[as.numeric(pathway_No)>=hitMin,]
write.xlsx(pathway,paste0(indir,"/Pathway.xlsx"))
subdir<-paste0(indir,"/PIC")
if (!dir.exists(subdir)){
  dir.create(subdir)
}
#图片下载
for (i in 1:length(nodeslink[4:No])){
  numHit<-nodeslink[4:No][i] %>% html_nodes("a") %>%  `[`(2) %>% html_text() %>% as.numeric()
  if (numHit>=hitMin){
    urlRaw<-nodeslink[4:No][i] %>% html_node("a") %>% html_attr("href")
    urlHead<-"https://www.kegg.jp"
    if (is_empty(grep(urlHead,urlRaw))){
      picurl<-paste0(urlHead,urlRaw)
    }
    else {picurl<-urlRaw}
    pic<-read_html(picurl)
    picDiv<-html_nodes(pic,"div.image-block") %>% html_nodes("img") %>% html_attr("src")
    if (is_empty(grep(urlHead,picDiv))){
      picture<-paste0(urlHead,picDiv)
    }
    else {picture<-picDiv}
    picName<-basename(picture)
    picName<-gsub("_[0-9]\\.[0-9]*","",picName)
    cat(paste0("------","download the picture of ",i,"/",nrow(pathway),"------\n"))
    download.file(picture,paste0(indir,"/PIC/",picName),mode = "wb",method = "wininet",cacheOK = T)
    Sys.sleep(1)
  }
}
