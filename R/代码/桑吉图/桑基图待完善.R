library(riverplot)
library(openxlsx)
library(pals)

data<-read.xlsx("代谢-临床指标-偏相关.xlsx",sheet=6)
data$X1[data$X1=="TG"|data$X1=="TC"]<-"Lipids Disease"
data$X1[data$X1=="SCr"|data$X1=="eGFR"|data$X1=="24h.Urine.protein"]<-"Renal"
data$X1[data$X1=="ALT"|data$X1=="AST"]<-"Liver"
data$X1[data$X1=="C3"|data$X1=="C4"|data$X1=="Anti-ds-DNA"|data$X1=="ANA"]<-"Immunologic"
data$X1[data$X1=="Anti-β2GPI"|data$X1=="ACL-IgG"|data$X1=="ACL-IgM"]<-"Anti-phospholipid"
data$X1[data$X1=="WNC"|data$X1=="PLT"]<-"Hematologic"
data$X1[data$X1=="CRP"]<-"Inflammation"
data$X1[data$X1=="ESR"|data$X1=="Indirect.Bilirubin"|data$X1=="LDH"]<-"Others"
nodes<-c("Lipids Disease","Renal","Liver","Immunologic","Anti-phospholipid","Hematologic","Inflammation","Others",c(unique(data$Classification)))


edges<-data.frame()
a<-1

for(i in unique(data$X1)){
  for(j in unique(data$Classification)){
    value<-length(rownames(data[data$X1==i,][data[data$X1==i,]$Classification==j,]))
    if(value>0){
      edges[a,1]<-i
      edges[a,2]<-j
      edges[a,3]<-value
      a<-a+1
    }
  }
}
colnames(edges)<-c("N1","N2","Value")
nodes<-data.frame(ID = unique(c(edges$N1, edges$N2)), stringsAsFactors = FALSE)
nodes$x<-c(rep(1,8),rep(2,10))
palette<-jet(18)
styles<-lapply(1:18, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})
names(styles)<-nodes$ID
# names(edges)[names(edges)=="TG"|names(edges)=="TC"]<-"Lipids Disease"
# names(edges)[names(edges)=="SCr"|names(edges)=="eGFR"|names(edges)=="24h.Urine.protein"]<-"Renal"
# names(edges)[names(edges)=="ALT"|names(edges)=="AST"]<-"Liver"
# names(edges)[names(edges)=="C3"|names(edges)=="C4"|names(edges)=="Anti-ds-DNA"|names(edges)=="ANA"]<-"Immunologic"
# names(edges)[names(edges)=="Anti-β2GPI"|names(edges)=="ACL-IgG"|names(edges)=="ACL-IgM"]<-"Anti-phospholipid"
# names(edges)[names(edges)=="WNC"|names(edges)=="PLT"]<-"Hematologic"
# names(edges)[names(edges)=="CRP"]<-"Inflammation"
# names(edges)[names(edges)=="ESR"|names(edges)=="Indirect.Bilirubin"|names(edges)=="LDH"]<-"Others"

# styles<-list("Lipids Disease"=list(col="#00007F"),Inflammation=list(col="#00E8FF"),Xenobiotics=list(col="#FFAC00"),
#              Renal=list(col="#0000BB"),Others=list(col="#25FFD9"),Carbohydrates=list(col="#FF7000"),
#              Liver=list(col="#0000F7"),"Amino acids"=list(col="#61FF9D"),Peptides=list(col="#FF3400"),
#              Immunologic=list(col="#0034FF"),Lipids=list(col="#9DFF61"),"Microbial metabolite"=list(col="#F70000"),
#              "Anti-phospholipid"=list(col="#0070FF"),Nucleotides=list(col="#D9FF25"),Energy=list(col="#BB0000"),
#              Hematologic=list(col="#00ACFF"),"Others Metabolites"=list(col="#FFE800"),"Cofactors and vitamins"=list(col="#7F0000"))
# r<-makeRiver(nodes,edges,node_xpos = c(rep(1,8),rep(2,10)),node_styles = list("Lipids Disease"=list(col="#00007F"),Inflammation=list(col="#00E8FF"),Xenobiotics=list(col="#FFAC00"),
                                                                           #     Renal=list(col="#0000BB"),Others=list(col="#25FFD9"),Carbohydrates=list(col="#FF7000"),
                                                                           #     Liver=list(col="#0000F7"),"Amino acids"=list(col="#61FF9D"),Peptides=list(col="#FF3400"),
                                                                           # Immunologic=list(col="#0034FF"),Lipids=list(col="#9DFF61"),"Microbial metabolite"=list(col="#F70000"),
                                                                           #     "Anti-phospholipid"=list(col="#0070FF"),Nucleotides=list(col="#D9FF25"),Energy=list(col="#BB0000"),
                                                                           #     Hematologic=list(col="#00ACFF"),"Others Metabolites"=list(col="#FFE800"),"Cofactors and vitamins"=list(col="#7F0000")))

rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
riverplot(rp)
plot(rp, plot_area = 0.95, yscale=0.06)

