library(ggplot2)
library(stringi)
library(stringr)
data <- read.csv('01.节律与差异节律/keeplogcpm.csv',row.names = 1)
clkgs <- c('Ciart','Clock','Cry1','Cry2','Dbp','Npas2','Nr1d1','Nr1d2','Per1','Per2','Per3')
# clkgs <- c('Cd73','Nadk','Nadsyn1','Nampt','Naprt','Nmnat3','Nnt','Nrk1','Paps','Pnp','Qprt','Sirts')
# clkgs <- c('Acadl','Acsl4','Cact','Cd36','Cpt1a','Cpt2','Crat','Fatp2')
clkdata <- data[clkgs,]
clkdata <- clkdata[!is.na(clkdata$Ctrl_ZT0_1),]
gpdata <- as.data.frame(as.table(as.matrix(clkdata)))
colnames(gpdata) <- c('symbol','sample','exp')
gpdata$time <- as.numeric(str_remove(str_split_i(gpdata$sample,'_',2),'ZT'))
gpdata$group <- str_split_i(gpdata$sample,'_',1)
# gpdata$time <- factor(gpdata$time,levels = unique(gpdata$time))
for(i in rownames(clkdata)){
  dat <- gpdata[gpdata$symbol==i,]
  # mean <- dat%>% group_by(group,time)%>%summarise(mean_val=mean(exp))
  p <- ggplot(data = dat,aes(x = time,y = exp,color=group))+
    # geom_point(data = dat,aes(x = time,y = exp,color=group,fill=group))+
    geom_point()+
    scale_color_manual(values = c('Ctrl'='black','DM'='red'))+
    # geom_line(data = mean,aes(x = time,y = mean_val,color=group))+
    geom_smooth(method = "loess",se = FALSE)+
    theme_bw()+
    theme(#panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      plot.title = element_text(color="black", size=20,hjust = 0.5),
      text=element_text(family="sans"),
      axis.ticks.x = element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.x = element_text(hjust = 1,size = 15,color = 'black'),
      axis.text.y = element_text(size = 15,color = 'black'),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      panel.border = element_rect(fill=NA,color="black",
                                  size=1, linetype="solid"))+
    labs(title=i,x='Time(hours)',y='logCPM')+
    scale_x_discrete(limits = c(0,4,8,12,16,20))
  ggsave(paste0(i,'.pdf'),p,width = 5,height = 5,dpi = 600)
  ggsave(paste0(i,'.png'),p,width = 5,height = 5,dpi = 600)
}

mean <- economics_long %>% 
  subset(variable %in% c("uempmed", "unemploy"))
mean[mean$date=='1967-07-01',]
