library(ropls)
library(ggplot2)
library(data.table)###多线程读文件
library(pals)
library(stringr)
library(psych)####相关性分析
library(openxlsx)
library(ggrepel)
options(stringsAsFactors = F,timeout = 200)

ModelScorePlot <- function(model, class.id, output,colanno,ci=0.95,ellipse=T,labelshow=F) {
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
      scale_color_manual(  # 颜色参数
        values = colanno
      ) +
      scale_fill_manual(
        values = colanno
      ) +
      geom_label_repel(
        data = df.p, 
        mapping = aes(x, y, label = rownames(df.p)), 
        color = 'black', size = 2, 
        label.padding = unit(0.15, 'lines'), 
        point.padding = unit(0.5, 'lines'), 
        min.segment.length = unit(0.1, "lines"), 
        segment.color = 'grey50', segment.size = 0.8, 
        max.overlaps = getOption("ggrepel.max.overlaps", default = 200),
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
        labs(x = paste0('PC1(',round(100*model@pcaVarVn['p1']/sum(model@pcaVarVn),2),'%)'),
             y = paste0('PC2(',round(100*model@pcaVarVn['p2']/sum(model@pcaVarVn),2),'%)'))
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
  # pdf(paste0(output, '/', model@typeC, ' score plot.pdf'))
  # p
  # dev.off()
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

ModelScorePlot <- function(model, class.id, output,colanno,ci=1,ellipse=F,labelshow=F) {
  # PCA/OPLS-DA模型得分图
  # 
  df.p <- as.data.frame(model$layout)
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
      scale_color_manual(  # 颜色参数
        values = colanno
      ) +
      scale_fill_manual(
        values = colanno
      ) +
      geom_label_repel(
        data = df.p, 
        mapping = aes(x, y, label = rownames(df.p)), 
        color = 'black', size = 2, 
        label.padding = unit(0.15, 'lines'), 
        point.padding = unit(0.5, 'lines'), 
        min.segment.length = unit(0.1, "lines"), 
        segment.color = 'grey50', segment.size = 0.8, 
        max.overlaps = getOption("ggrepel.max.overlaps", default = 200),
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
        labs(x = paste0('PC1(',round(100*model@pcaVarVn['p1']/sum(model@pcaVarVn),2),'%)'),
             y = paste0('PC2(',round(100*model@pcaVarVn['p2']/sum(model@pcaVarVn),2),'%)'))
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
      
      labs(x = "UMAP_1",
            y = "UMAP_2")
  }
  # pdf(paste0(output, '/', model@typeC, ' score plot.pdf'))
  # p
  # dev.off()
  ggsave(  # 输出位置
    paste0(output, '/', ' umap.jpg'), p,
    width = 9.6, height = 6, units = 'in', dpi = 600
  )
  ggsave(
    paste0(output, '/', ' umap.pdf'), p,
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
  write.xlsx(df.data,file = paste0(output, '/', model@typeC, '_3D_data.xlsx'),row.names=TRUE)
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



raw<-read.xlsx('Supplementary table S2.xlsx',sheet = 3,rowNames = T)
# data<-raw[,c(grep('Control|SLE',colnames(raw)))]
# rownames(data)<-1:nrow(data)
# class.id=factor(sapply(colnames(data),function(i){paste0(str_extract_all(i,'\\D')[[1]],collapse = '')}),
#                 levels = unique(sapply(colnames(data),function(i){paste0(str_extract_all(i,'\\D')[[1]],collapse = '')})))
data<-raw
class.id<-data$Drug
names(class.id)<-rownames(data)
data<-data[,-1]
class.id<-factor(class.id,levels = unique(class.id))

colanno<-jet(length(unique(class.id))+2)[1:length(unique(class.id))+1]
names(colanno)<-unique(class.id)





pca <- opls(
  x = data, 
  log10L = T, scaleC = 'standard',#predI = 2,
  crossvalI = min(ncol(data), 7)
)
# pca<-umap(t(data))

ModelScorePlot(  # 2D得分图
  model = pca, class.id = class.id,colanno = colanno,
  output = "./")
ModelScorePlot3DNew(model = pca, class.id = class.id,colanno = colanno,
                    output = "D:/马驰宇/桌面/何院专辑/研究院项目/戴主任项目重分析/7月份的PCA")



pls.da <- opls(
  y = as.character(class.id),permI = 200, 
  x = t(data), log10L = T, scaleC = 'standard', 
  crossvalI = min(ncol(data), 7), 
  printL = F, plotL = F
)

ModelScorePlot(  # 2D得分图
  model = pls.da, class.id = class.id,colanno = colanno,
  output = 'all')
ModelScorePlot3DNew(model = pls.da, class.id = class.id,colanno = colanno,
                    output = 'all')

					
opls.da<-opls(
  y = as.character(class.id),permI = 200,
  x = t(data), predI = 1, ortho = 1, log10L = T, scaleC = 'standard',
  crossvalI = min(ncol(data), 7),
  printL = F, plotL = F
)  

ModelScorePlot(  # 2D得分图
  model = opls.da, class.id = class.id,colanno = colanno,
  output = 'all')
ModelScorePlot3DNew(model = opls.da, class.id = class.id,colanno = colanno,
                    output = 'all')