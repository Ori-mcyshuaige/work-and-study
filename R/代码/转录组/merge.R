

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("sva")


#???ð?
library(limma)
library(sva)
files=c("batch_GSE98793.csv", "batch_GSE76826.csv",'batch_GSE32280.csv')                  #?????ļ?????
path='./01.差异分析/'
#??ȡ????????
geneList=list()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.csv(paste0(path,inputFile), header=T, check.names=F)
    rt=rt[-nrow(rt),]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

#???ݺϲ?
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    #??ȡ?????ļ????????????ļ?????????
    rt=read.csv(paste0(path,inputFile), header=T, check.names=F)
    colnames(rt) <- paste0(colnames(rt),'_',rt[unlist(nrow(rt)),])
    rt=rt[-nrow(rt),]
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))
    #??TCGAɾ????????Ʒ
    if(header[1] == "TCGA"){
		group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
		group=sapply(strsplit(group,""), "[", 1)
		rt=rt[,group==0]
		rt=t(rt)
		row.names(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(rt))
		rt=avereps(rt)
		rt=t(rt)
    }
    #????ֵ????????ȡlog2
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    if(header[1] != "TCGA"){
    	rt=normalizeBetweenArrays(rt)
    }
    #???ݺϲ?
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#?????ݽ??н????????????????Ľ???
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(outTab,group=sapply(1:184,function(x){str_split(colnames(outTab),'_')[[x]][4]}))
outTab1=rbind(geneNames=colnames(outTab), outTab)

for(i in files){
  var=str_split(i,'\\.')[[1]][1]
  data=outTab[,grep(var,colnames(outTab1))]
  assign(var,data)
  write.csv(data,paste0(var,'.csv'))
}
write.table(outTab1, file="merge.txt", sep="\t", quote=F, col.names=F)



