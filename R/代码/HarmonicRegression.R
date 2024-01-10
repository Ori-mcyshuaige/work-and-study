library(HarmonicRegression)

raw <- read.csv('01.节律与差异节律/keeplogcpm.csv',row.names = 1)

ctrl <- raw[,grep('Ctrl',colnames(raw))]
ctrl <- ctrl[,c(grep('ZT0',colnames(ctrl)),grep('ZT4',colnames(ctrl)),grep('ZT8',colnames(ctrl)),grep('ZT12',colnames(ctrl)),
                      grep('ZT16',colnames(ctrl)),grep('ZT20',colnames(ctrl)))]
dm <- raw[,grep('DM',colnames(raw))]
dm <- dm[,c(grep('ZT0',colnames(dm)),grep('ZT4',colnames(dm)),grep('ZT8',colnames(dm)),grep('ZT12',colnames(dm)),
            grep('ZT16',colnames(dm)),grep('ZT20',colnames(dm)))]

ctrl_rs <- harmonic.regression(t(as.matrix(ctrl)), 
                               inputtime = c(0,0,0,0,4,4,4,4,8,8,8,12,12,12,12,16,16,16,16,20,20,20,20), 
                               Tau = 24, normalize = TRUE,norm.pol = FALSE, 
                               norm.pol.degree = 1, trend.eliminate = FALSE,trend.degree = 1)
# normts <- ctrl_rs$normts
# fitvals <- ctrl_rs$fit.vals
# normfitvals <- ctrl_rs$norm.fit.vals
ctrl_pars <- ctrl_rs$pars
ctrl_pars$p <- ctrl_rs$pvals[rownames(ctrl_pars)]
ctrl_pars$q <- ctrl_rs$qvals[rownames(ctrl_pars)]
write.csv(ctrl_pars,'ctrl_amp.csv')

dm_rs <- harmonic.regression(t(as.matrix(dm)), 
                               inputtime = c(0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12,16,16,16,16,20,20,20,20), 
                               Tau = 24, normalize = TRUE,norm.pol = FALSE, 
                               norm.pol.degree = 1, trend.eliminate = FALSE,trend.degree = 1)
# normts <- ctrl_rs$normts
# fitvals <- ctrl_rs$fit.vals
# normfitvals <- ctrl_rs$norm.fit.vals
dm_pars <- dm_rs$pars
dm_pars$p <- dm_rs$pvals[rownames(dm_pars)]
dm_pars$q <- dm_rs$qvals[rownames(dm_pars)]
write.csv(dm_pars,'dm_amp.csv')

allamp <- data.frame(ctrlamp=ctrl_pars$amp,dmamp=dm_pars[rownames(ctrl_pars),]$amp)
rownames(allamp) <- rownames(ctrl_pars)
allamp$group <- 'noncircadian'
circadian <- read.csv('ctrl_rhygene.csv')
cirgene <- circadian$X
allamp[cirgene,]$group <- 'circadian'
allamp$log2Arel <- log2(allamp$dmamp/allamp$ctrlamp)
write.csv(allamp,'allampstat.csv')

data <- data.frame(x=rep(c('Ctrl','Dm'),each=length(allamp$ctrlamp)),y=c(allamp$ctrlamp,allamp$dmamp))

p<-ggplot(data,aes(x,y),size=0.4)+
  stat_boxplot(geom = "errorbar",size=0.4,width = 0.2,aes(color = x))+
  geom_boxplot(aes(color = x),position = position_dodge(width = 0.2),notch=TRUE,notchwidth=0.8,size = 0.4,fill = "white",width = 0.4,outlier.alpha = 0)+
  # geom_jitter(aes(color = Group),width =0.175,shape = 16,size=2)+
  geom_point(size=1,aes(color = x,fill = x))+
  # geom_hline(yintercept = 0,color = "grey",linetype = "dashed")+
  # geom_line(aes(group=yname) ,size=0.5,color="grey50")+
  scale_color_manual(values = c(Ctrl="black",Dm="red"))+
  # scale_fill_manual(values = c(Control="#9B30FF",SLE="#FF8247"))+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10,colour = "black"),
    axis.title = element_text(size = 10,colour = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 1,colour = "black")
  )+
  #   compare_means(len ~ supp, data = ToothGrowth, 
  #                 group.by = "dose", paired = TRUE)
  # # Box plot facetted by "dose"
  #   ggpaired(ToothGrowth, x = "supp", y = "len",
  #               color = "supp", palette = "jama", 
  #               line.color = "gray", line.size = 0.4,
  #               facet.by = "dose", short.panel.labs = FALSE)+
  stat_compare_means(comparisons = list(c(data$x[1],data$x[nrow(data)])), method = "wilcox.test", paired = F,label = "p.signif")+
  # stat_compare_means()+
  labs(x=NULL,y="Relative amplitude")
# coord_cartesian(ylim = c(0,1.1*max(datas[,i])))
# scale_y_continuous(expand = c(0, 0))
ggsave("rel amp.png", p, width = 4.8, height = 4.8,units = "in",dpi = 600)
ggsave("rel amp.pdf", p, width = 4.8, height = 4.8,dpi = 600)


