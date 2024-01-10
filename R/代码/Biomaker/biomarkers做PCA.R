library(ropls)
library(ggplot2)
library(openxlsx)
library(pals)
library(stringr)

data_rf<-read.xlsx('data_rf.xlsx',rowNames = T)
data_yz<-read.xlsx('data_yz.xlsx',rowNames = T)

data_all<-rbind(data_rf,data_yz[,colnames(data_rf)])

ModelScorePlot <- function(model, class.id, output,colanno,ci=0.95,ellipse=F,labelshow=F) {
  # PCA/OPLS-DA模型得分图
  # 
  if(model@typeC == 'PCA' | model@typeC == 'PLS-DA') {
    df.p <- as.data.frame(model@scoreMN[, 1:2])
  } else if(model@typeC == 'OPLS-DA') {
    df.p <- cbind(
      as.data.frame(model@scoreMN[, 1]), 
      as.data.frame(model@orthoScoreMN[,1])
    )
  }
  colnames(df.p) <- c('x', 'y')
  # write.xlsx(df.p,file = paste0(output, '/', model@typeC, 'data.xlsx'),row.names=TRUE)
  n <- nrow(df.p)
  hfn <- 2*(n-1)*(n^2-1)/(n^2*(n-2))*qf(ci, 2, (n-2))
  rv <- seq(0, 2*pi, length.out = 100)
  # 置信区间数据
  df.ell <- data.frame(  
    x = sqrt(var(df.p$x)*hfn)*cos(rv), 
    y = sqrt(var(df.p$y)*hfn)*sin(rv)
  )
  p <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_polygon(data = df.ell, aes(x, y), color = 'black', fill = NA)
  if(ellipse) {
    p <- p + stat_ellipse(
      data = df.p, geom = 'polygon', level = ci, 
      aes(x, y, fill = class.id, color = class.id), alpha = I(0.1)
    ) +
      stat_ellipse(
        data = df.p, geom = 'blank', level = ci, 
        aes(-x, -y, fill = class.id)
      )
  }
  if(labelshow) {
    p <- p + geom_point(
      data = df.p, size = 6, aes(x, y, shape = class.id, color = class.id)
    ) +
      geom_blank(
        data = df.p, aes(-x, -y, shape = class.id, color = class.id)
      ) +
      scale_shape_manual(  # 形状参数
        values = rep_len(
          #c("\u2605","\u25C4","\u25BC","\u25B2"), 
          c(16, 15, 17, 18),
          length(levels(class.id))
        )[sort(unique(as.numeric(class.id)))]
      ) +
      # 颜色参数
      scale_color_manual( 
        values = colanno
      ) +
      scale_fill_manual(
        values = colanno
      ) +
      geom_label_repel(
        data = df.p, 
        mapping = aes(x, y, label = rownames(df.p)), 
        color = 'black', size = 3, 
        label.padding = unit(0.2, 'lines'), 
        point.padding = unit(0.5, 'lines'), 
        min.segment.length = unit(0.1, "lines"), 
        segment.color = 'grey50', segment.size = 1, 
        show.legend = F
      ) +
      theme_bw() +
      theme(
        legend.title = element_blank(), 
        legend.key = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 1)
        # axis.text = element_text(size = 15, color = "black",
        #                          face = "plain", vjust = 0.5, hjust = 0.5),
        # axis.title = element_text(size = 15, color = "black", 
        #                           face = "plain", vjust = 0.5, hjust = 0.5)
      ) +
      
      if(model@typeC == 'PCA') {
        labs(x = 'PC[1]', y = 'PC[2]')
      } else if(model@typeC == 'OPLS-DA') {
        labs(x = 't[1]P', y = 't[1]O')
      } else if(model@typeC == 'PLS-DA') {
        labs(x = 't[1]', y = 't[2]')
      }
  }
  if(labelshow == F) {
    p <- p + geom_point(
      data = df.p, size = 5, aes(x, y, shape = class.id, color = class.id)
    ) +
      geom_blank(
        data = df.p, aes(-x, -y, shape = class.id, color = class.id)
      ) +
      scale_shape_manual(  # 形状参数
        values = rep_len(
          #c("\u2605","\u25C4","\u25BC","\u25B2"), 
          c(16, 15, 17, 18),
          length(levels(class.id))
        )[sort(unique(as.numeric(class.id)))]
      ) +
      # 颜色参数
      scale_color_manual(  
        values = colanno
      ) +
      scale_fill_manual(
        values = colanno
      ) +
      theme_bw() +
      theme(
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 1)
        # axis.text = element_text(size = 15, color = "black",
        #                          face = "plain", vjust = 0.5, hjust = 0.5),
        # axis.title = element_text(size = 15, color = "black", 
        #                           face = "plain", vjust = 0.5, hjust = 0.5)
      ) +
      if(model@typeC == 'PCA') {
        labs(x = paste0('PC1(',round(100*model@pcaVarVn['p1']/sum(model@pcaVarVn),2),'%)'),
             y = paste0('PC2(',round(100*model@pcaVarVn['p2']/sum(model@pcaVarVn),2),'%)'))
      } else if(model@typeC == 'OPLS-DA') {
        labs(x = 't[1]P', y = 't[1]O')
      } else if(model@typeC == 'PLS-DA') {
        labs(x = 't[1]', y = 't[2]')
      }
  }
  ggsave(  # 输出位置
    paste0(output, '/', model@typeC, ' score plot.jpg'), p, 
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', model@typeC, ' score plot.pdf'), p,
    device = cairo_pdf,
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
}


ModelScorePlot3DNew <- function(model, class.id, output,colanno) {
  # 新的3D PCA得分图~~撒花~
  # class.id参数必须是因子，不然会出现图例的排序问题
  # 
  require(plotly)
  require(htmlwidgets)
  df.data <- as.data.frame(model@scoreMN[, 1:3])
  df.data <- cbind(class.id, df.data)
  colnames(df.data) <- c("class.id", "PC1", "PC2", "PC3")
  # write.xlsx(df.data,file = paste0(output, '/', model@typeC, '_3D_data.xlsx'),row.names=TRUE)
  axis.x <- list(
    range = c(- max(abs(df.data$PC1)), max(abs(df.data$PC1)))
  )
  axis.y <- list(
    range = c(- max(abs(df.data$PC2)), max(abs(df.data$PC2)))
  )
  axis.z <- list(
    range = c(- max(abs(df.data$PC3)), max(abs(df.data$PC3)))
  )
  p <- plot_ly(
    df.data, x = ~PC1, y = ~PC2, z = ~PC3, 
    type = "scatter3d", mode = "markers", 
    color = ~class.id, 
    colors = colanno, 
    symbol = ~class.id, 
    symbols = rep_len(
      #c("\u2605","\u25C4","\u25BC","\u25B2"),
      c(16, 15, 17, 18),
      length.out = length(levels(class.id))
    ), 
    text = ~row.names(df.data), 
    marker = list(size = 8, opacity = 0.8)
  ) %>% 
    layout(
      legend = list(
        bgcolor = "#FFFFFF",
        bordercolor = "#000000"
        
      ), 
      scene = list(
        xaxis = axis.x, yaxis = axis.y, zaxis = axis.z
      )
    )
  old.wd <- getwd()
  setwd(output)
  on.exit(setwd(old.wd))
  htmlwidgets::saveWidget(p, file = paste0(model@typeC," score plot 3D.html"))
}



for (file in c(Sys.glob('*/*/biomarkers.xlsx'),Sys.glob('*/biomarker.xlsx'))) {
  biomarkers<-read.xlsx(file,rowNames = T)
  outdir0<-paste0(strsplit(file,'/')[[1]][1],'/','7 PCA')
  if(!dir.exists(outdir0)){dir.create(outdir0)}
  
  
  
  # analysis
  data<-data[,colnames(biomarkers)[-ncol(biomarkers)]]
  groups<-sapply(rownames(data), function(i){strsplit(i,'_')[[1]][1]})
  groups[groups=='C']<-'Control'
  class.id=factor(groups,levels = unique(groups))
  colanno<-jet(length(unique(class.id))+2)[1:length(unique(class.id))+1]
  names(colanno)<-unique(class.id)
  outdir<-paste0(outdir0,'/analysis')
  if(!dir.exists(outdir)){dir.create(outdir)}
  pca <- opls(
    x = data, 
    log10L = T, scaleC = 'center', 
    crossvalI = min(ncol(data), 7), 
    printL = F, plotL = F
  )
  if(pca@summaryDF$pre<3){
    pca <- opls(
      x = data, predI = 3,
      log10L = T, scaleC = 'center', 
      crossvalI = min(ncol(data), 7), 
      printL = F, plotL = F
    ) 
  }
  df.s<- as.data.frame(pca@scoreMN)
  df.l<- as.data.frame(pca@loadingMN)
  write.xlsx(df.s,file = paste0(outdir, '/', pca@typeC, '_scoreMN.xlsx'),row.names=TRUE)
  write.xlsx(df.l,file = paste0(outdir, '/', pca@typeC, '_loadingMN.xlsx'),row.names=TRUE)
  write.csv(pca@pcaVarVn,file = paste0(outdir, '/', pca@typeC, '_value.csv'),row.names=TRUE)
  ModelScorePlot(  # 2D得分图
    model = pca, class.id = class.id,colanno = colanno,
    output = outdir)
  ModelScorePlot3DNew(model = pca, class.id = class.id,colanno = colanno,
                      output = outdir)
  
  
  

  
    # all
  data<-data_all[,colnames(biomarkers)[-ncol(biomarkers)]]
  groups<-sapply(rownames(data), function(i){strsplit(i,'_')[[1]][1]})
  groups[groups=='C']<-'Control'
  class.id=factor(groups,levels = unique(groups))
  colanno<-jet(length(unique(class.id))+2)[1:length(unique(class.id))+1]
  names(colanno)<-unique(class.id)
  outdir<-paste0(outdir0,'/all')
  if(!dir.exists(outdir)){dir.create(outdir)}
  pca <- opls(
    x = data, 
    log10L = T, scaleC = 'center', 
    crossvalI = min(ncol(data), 7), 
    printL = F, plotL = F
  )
  if(pca@summaryDF$pre<3){
    pca <- opls(
      x = data, predI = 3,
      log10L = T, scaleC = 'center', 
      crossvalI = min(ncol(data), 7), 
      printL = F, plotL = F
    ) 
  }
  df.s<- as.data.frame(pca@scoreMN)
  df.l<- as.data.frame(pca@loadingMN)
  write.xlsx(df.s,file = paste0(outdir, '/', pca@typeC, '_scoreMN.xlsx'),row.names=TRUE)
  write.xlsx(df.l,file = paste0(outdir, '/', pca@typeC, '_loadingMN.xlsx'),row.names=TRUE)
  write.csv(pca@pcaVarVn,file = paste0(outdir, '/', pca@typeC, '_value.csv'),row.names=TRUE)
  ModelScorePlot(  # 2D得分图
    model = pca, class.id = class.id,colanno = colanno,
    output = outdir)
  ModelScorePlot3DNew(model = pca, class.id = class.id,colanno = colanno,
                      output = outdir)
  
  
}
