library(openxlsx)
library(ggplot2)

###重要补充，调整标签位置
#output<-output[,label_y:=cumsum(Number)-Number/2,by=list(Tissue)]
output$label<-output$Number/2
output$label[output$color=="down"]<-output$label[output$color=="down"]+output$Number[output$color=="up"]

setwd("D:/马驰宇/181/戴主任/工作/生信流程/项目/分析/已交付修回-辜美佳 BQ-GMJ20200619-GZ1-HFX-RXQ/修回11.2-BQ-GMJ20200619-GZ1-HFX-RXQ")
neg<-read.xlsx("柱形图neg.xlsx",sheet = 1,colNames = T,rowNames = F,check.names = F,sep.names = " ")
pos<-read.xlsx("柱形图pos.xlsx",sheet = 1,colNames = T,rowNames = F,check.names = F,sep.names = " ")
colnames(neg)<-"mcy"
colnames(pos)<-"mcy"
dataneg<-rbind(neg[,c(1,4)],neg[,c(1,2)])
dataneg<-rbind(dataneg,neg[,c(1,5)])
dataneg<-rbind(dataneg,neg[,c(1,3)])
dataneg$group<-"Normal"
dataneg[,"group"][c(22:42)]<-"Ordinary"
dataneg[,"group"][c(43:63)]<-"Critical"
dataneg[,"group"][c(64:84)]<-"Severe"
datapos<-rbind(pos[,c(1,4)],pos[,c(1,2)])
datapos<-rbind(datapos,pos[,c(1,5)])
datapos<-rbind(datapos,pos[,c(1,3)])
datapos$group<-"Normal"
datapos[,"group"][c(57:112)]<-"Ordinary"
datapos[,"group"][c(113:168)]<-"Critical"
datapos[,"group"][c(169:224)]<-"Severe"
colnames(dataneg)<-c("MS2 Name","Mean","Group")
colnames(datapos)<-c("MS2 Name","Mean","Group")
datapos$Group<-factor(datapos$Group,levels = c("Normal","Ordinary","Critical","Severe"))
dataneg$Group<-factor(dataneg$Group,levels = c("Normal","Ordinary","Critical","Severe"))

p<-ggplot(data = output,aes(x = Number, y = Tissue,color=color,fill=color))+
  geom_bar(stat = "identity",position = "stack",width=0.5)+  ###position = "dodge" 列  “stack”堆叠   "fill" 百分比
  scale_color_manual(values = c("up"="#FF3030","down"="#4876FF"))+
  geom_text(aes(x=label,y=Tissue,label=Number),size=7,color="black",position='identity')+
  scale_fill_manual(values = c("up"="#FF3030","down"="#4876FF"))+
  theme_bw()+
  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black",size = 1),
  legend.text = element_text(size = 16,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
  legend.title = element_blank(),
  # legend.position = "none",
  # legend.key.size = unit(0.6,"cm"),
  axis.title = element_text(size = 14,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
  axis.text.x = element_text(size = 14,color = "black",family = "sans",face = "plain",angle = 0, hjust = 0.5,vjust = 0.5),
  axis.text.y = element_text(size = 14,color = "black",family = "sans",face = "plain", hjust = 1,vjust = 0)
)+
  scale_y_continuous(limits = c(0,58),expand = c(0,0))+
  # coord_flip()+
  labs(x = "Pro Number", y = "Tissue")

ggsave("上下调柱形图竖.png",p,dpi = 600,width = 10,height = 20)
ggsave("上下调柱形图竖.pdf",p,width = 10,height = 20)


p<-ggplot(data = data,aes(x = X.NAME., y = X,color=Impact,fill=Impact))+
  # geom_blank(aes(x=-X.NAME.,y=X))+
  geom_bar(stat = "identity",position = "dodge",width = 0.5)+  ###position = "dodge" 列  “stack”堆叠   "fill" 百分比
  scale_color_gradient(
    name = "Impact", 
    low = '#8B7E66', high = '#A52A2A', 
    guide = guide_colorbar(
      title.position = 'top', title.hjust = 0.5, 
      direction = 'horizontal'
    ))+
      scale_fill_gradient(
        name = "Impact", 
        low = '#8B7E66', high = '#A52A2A', 
        guide = guide_colorbar(
          title.position = 'top', title.hjust = 0.5, 
          direction = 'horizontal'
        ))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",size = 1),
    legend.text = element_text(size = 12,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_text(size = 12,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.key.size = unit(0.6,"cm"),
    axis.title = element_text(size = 12,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text.x = element_text(size = 12,color = "black",family = "sans",face = "plain", hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 12,color = "black",family = "sans",face = "plain", hjust = 1,vjust = 0),
    legend.position = c(0.85,0.1),
  )+
  # coord_flip()+
  labs(x = "-Log10(P)", y = "Pathway")

ggsave("Pathway fc.png",p,dpi = 600,width = 9,height = 7)
ggsave("Pathway fc.pdf",p,width = 9,height = 7)
  
data<-data[order(data$X.NAME.),]


p<-ggplot(data = data,aes(x = x, y = Term,color=X2,fill=X2))+
  # geom_blank(aes(x=-Term,y=x))+
  geom_bar(stat = "identity",position = "dodge",width = 0.5)+  ###position = "dodge" 列  “stack”堆叠   "fill" 百分比
  scale_color_manual(
    values = color
  )+
  scale_fill_manual(
    values = color
  )+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black",size = 1),
    legend.text = element_text(size = 10,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.title = element_text(size = 10,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0),
    legend.key.size = unit(0.6,"cm"),
    axis.title = element_text(size = 10,color = "black",family = "sans",face = "plain",vjust = 0.5, hjust = 0.5),
    axis.text.x = element_text(size = 10,color = "black",family = "sans",face = "plain", hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 10,color = "black",family = "sans",face = "plain", hjust = 1,vjust = 0),
    legend.position = c(0.7,0.1)
  )+
  # coord_flip()+
  labs(x = "-Log10(q)", y = "Pathway")
