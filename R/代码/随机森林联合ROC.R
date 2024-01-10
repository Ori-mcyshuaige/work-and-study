pdf('testROC plot.pdf',width = 7,height = 7)
plot.roc(testmet,col='orange',grid=F,print.auc=F,legacy.axes=T)
plot.roc(testpro,col='skyblue',grid=F,print.auc=F,legacy.axes=T,add = T)
plot.roc(testpromet,col='purple',grid=F,print.auc=F,legacy.axes=T,add = T)
legendtest<-c('DATA: AUC(95%CI) %',
              paste0('compound: ',sprintf("%.2f (%.2f-%.2f) %%", testmet$auc, testmet$ci[1],testmet$ci[3])),
              paste0('protein: ',sprintf("%.2f (%.2f-%.2f) %%", testpro$auc, testpro$ci[1],testpro$ci[3])),
              paste0('all: ',sprintf("%.2f (%.2f-%.2f) %%", testpromet$auc, testpromet$ci[1],testpromet$ci[3]))
)
legendcol<-c("black",'orange','skyblue',"purple")
# if (length(unique(unlist(sapply(rownames(liqundataraw),function(i){substring(paste0(str_extract_all(i,'\\D')[[1]],collapse = ''),1,2)}))))>1){
#   plot.roc(liqundataroc,print.thres=F,print.auc=F,col='green',reuse.auc=F,add = T)
#   legendtest<-c(legendtest,paste0('newdata: ',sprintf("%.2f (%.2f-%.2f) %%", liqundataroc$auc, liqundataroc$ci[1],liqundataroc$ci[3])))
#   legendcol<-c(legendcol,'green')
# }
legend('bottomright',
       legendtest,col=legendcol,text.col = legendcol,bty = 'n')


dev.off()
