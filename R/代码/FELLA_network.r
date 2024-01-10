network1 <- function() {
set.seed(123)
wb <- createWorkbook(creator = 'Biotree')
modifyBaseFont(wb, fontSize = 8, fontName = 'Arial')
hs <- createStyle(
  fontSize = 10, fontColour = 'white', fgFill = '#008CCE',
  valign = 'center', halign = 'center'
)
nlimit <- 200
# vertex.label.cex <- .8
for (i in 1:length(l.a$subDir)) {
# print(i)
# Filter overview pathways
tmpdir <- paste0("/home/share/networkDb/",l.a$keggOrg)

# Mke sure the database does not exist from a former vignette build
# Otherwise the vignette will rise an error
# because FELLA will not overwrite an existing database
if (!dir.exists(tmpdir) | (dir.exists(tmpdir) & file.info(tmpdir)$size==0)){
  graph <- buildGraphFromKEGGREST(
    organism = l.a$keggOrg,
    filter.path = c("01100", "01200", "01210", "01212", "01230")
  )
  unlink(tmpdir, recursive = TRUE)
  buildDataFromGraph(
    keggdata.graph = graph,
    databaseDir = tmpdir,
    internalDir = FALSE,
    matrices = "diffusion",
    normality = "diffusion",
    niter = 50
  )
}

fella.data <- loadKEGGdata(
databaseDir = tmpdir,
internalDir = FALSE,
loadMatrix = "diffusion"
)

compounds.epithelial <- read.table(paste0(l.a$out.temp, '/KEGG.in', 
                                   ifelse(is.na(pol), '', paste0(pol, '/')), 
                                   l.a$subDir[i], '.txt'))[,1]

analysis.epithelial <- defineCompounds(
  compounds = compounds.epithelial,
  data = fella.data
)
getInput(analysis.epithelial)
getExcluded(analysis.epithelial)
# analysis.epithelial
analysis.epithelial <- runDiffusion(
  object = analysis.epithelial,
  data = fella.data,
  approx = "normality"
)
g <- generateResultsGraph(object = analysis.epithelial,data = fella.data,nlimit = nlimit)
vertice <- get.vertex.attribute(g)
verticeListName <- names(vertice)
verticeAttr <- lapply(1:length(verticeListName), function(n){
  Attr <- eval(str2expression(paste0('V(g)$',verticeListName[n])))
  Attr <-as.data.frame(t(t(Attr)))
})
nodeInfo <- Reduce(cbind.data.frame,verticeAttr)
names(nodeInfo) <- verticeListName
linkInfo <- as.data.frame(get.edgelist(g))
names(linkInfo) <- c("from","to")
nodeInfo$info <- NA
nodeInfo$info[nodeInfo$com==1] <-"Pathway"
nodeInfo$info[nodeInfo$com==2] <-"Module"
nodeInfo$info[nodeInfo$com==3] <-"Enzyme"
nodeInfo$info[nodeInfo$com==4] <-"Reaction"
nodeInfo$info[nodeInfo$com==5 & nodeInfo$input==FALSE] <-"Compound"
nodeInfo$info[nodeInfo$com==5 & nodeInfo$input==TRUE] <-"Input compound"

network_plot<-function(path,node,link,lyt="fruchtermanreingold",edgeColor="gray50"){
  node<-as.data.frame(node)
  link<-as.data.frame(link)
  nodeNumber<-1:nrow(node)
  names(nodeNumber)<-node$name
  
  nodePoint<-matrix(ncol=2,nrow = nrow(link))
  for (i in 1:nrow(link)){
    nodePoint[i,1]<-nodeNumber[link[i,1]]
    nodePoint[i,2]<-nodeNumber[link[i,2]]
  }
  
  set.seed(123)
  net<-network(nrow(node),directed = FALSE,vertex.attrnames = node$name)
  m<-as.matrix.network.adjacency(net)
  plotcord <- switch(tolower(lyt),
                     "fruchtermanreingold" = data.frame(gplot.layout.fruchtermanreingold(m, NULL)),
                     "adj" = data.frame(gplot.layout.adj(m, NULL)),
                     "circle" = data.frame(gplot.layout.circle(m, NULL)),
                     "eigen" = data.frame(gplot.layout.eigen(m, NULL)),
                     "geodist" = data.frame(gplot.layout.geodist(m, NULL)),
                     "hall" = data.frame(gplot.layout.hall(m, NULL)),
                     "kamadakawai" = data.frame(gplot.layout.kamadakawai(m, NULL)),
                     "mds" = data.frame(gplot.layout.mds(m, NULL)),
                     "princoord" = data.frame(gplot.layout.princoord(m, NULL)),
                     "random" = data.frame(gplot.layout.random(m, NULL)),
                     "rmds" = data.frame(gplot.layout.rmds(m, NULL)),
                     "segeo" = data.frame(gplot.layout.segeo(m, NULL)),
                     "seham" = data.frame(gplot.layout.seham(m, NULL)),
                     "spring" = data.frame(gplot.layout.spring(m, NULL)),
                     "springrepulse" = data.frame(gplot.layout.springrepulse(m, NULL)),
                     "target" = data.frame(gplot.layout.target(m, NULL))
  )
  colnames(plotcord)<-c("X1","X2")
  plotcord$name <- node$name
  plotcord$info <- node$info
  plotcord$label <- node$label
  
  edglist<-nodePoint
  edges<-data.frame(plotcord[edglist[,1],],plotcord[edglist[,2],])
  names(edges)<-c("X1","Y1","name1","info1","label1","X2","Y2","name2","info2","label2")
  labelTypes<-unique(nodeInfo$info)
  shapeValues<-rep(19,length(labelTypes))
  names(shapeValues)<-labelTypes
  shapeValues[grepl("^Input compound",labelTypes)]<-22
  sizeValues<-rep(8,length(labelTypes))
  names(sizeValues)<-labelTypes
  sizeValues[grepl("Pathway",labelTypes)]<-16
  sizeValues[grepl("Module",labelTypes)]<-12
  colValues<-rep("#548B54",length(labelTypes))
  names(colValues)<-labelTypes
  colValues[grepl("Pathway",labelTypes)]<-"#CD0000"
  colValues[grepl("Module",labelTypes)]<-"#CD96CD"
  colValues[grepl("Enzyme",labelTypes)]<-"#FFA200"
  colValues[grepl("Reaction",labelTypes)]<-"#8DB6CD"
  pnet <- ggplot(data = plotcord, aes(X1, X2, label = name,shape=info,color=info,fill=info,size=info)) + 
      geom_segment(data = edges, aes(x = X1, y = Y1, xend = X2, yend = Y2),size = 0.4, inherit.aes = FALSE,color=edgeColor) +
      geom_point(aes(shape=info,color=info,fill=info,size=info)) +
      scale_shape_manual(values = shapeValues) +
      scale_size_manual(values = sizeValues) +
      scale_color_manual(values = colValues) +
      scale_fill_manual(values = colValues) +
      labs(x = "", y = "",size = "",shape = "",color = "",fill = "") + 
      geom_text_repel(
        aes(label = name),
        size = 8,
        color = "black",family="serif"
      ) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 20,color = "black",family = "serif"),
        legend.text= element_text(face = "bold", color = "black", size = 18,family = "serif"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks  = element_blank(),
        legend.position = "right"
      )
  if (nrow(linkInfo)>=100 & nrow(linkInfo)<250){
    pheight<-pwidth<-15
  } else {pheight<-pwidth<-20}
  ggsave(paste0(path,"/net.png"),plot = pnet,units = "in",width = pwidth,height = pheight,dpi = 600)
  ggsave(paste0(path,"/net.pdf"),plot = pnet,units = "in",width = pwidth,height = pheight,dpi = 600)
}
network_plot(path =paste0(l.a$output, '/results/',
             ifelse(is.na(pol), '', paste0(pol, '/')),
             'Network Analysis/',
             l.a$subDir[i]),node = nodeInfo,link = linkInfo)
# pdf(paste0(
#   l.a$output, '/results/',
#   ifelse(is.na(pol), '', paste0(pol, '/')),
#   'Network Analysis/',
#   l.a$subDir[i], '/network.pdf'),
  # width = 10.6, height = 9.6, pointsize = 5
# )
# 
# netplot <- plot(
#   analysis.epithelial,
#   method = "diffusion",
#   data = fella.data,
#   nlimit = nlimit,
#   # LabelLengthAtPlot = 45,
#   NamesAsLabels = FALSE,
#   # font = 4,
#   plotLegend = T,
#   # graph.layout = layout_randomly(),
#   vertex.label.cex = .5,vertex.label.color="black"
#   # layout = TRUE,
#   # vertex.label.cex = vertex.label.cex
# )
# 
# dev.off()
# png(paste0(
#   l.a$output, '/results/', 
#   ifelse(is.na(pol), '', paste0(pol, '/')), 
#   'Network Analysis/', 
#   l.a$subDir[i], '/network.png'),
#   width = 6.6, height = 4.2, units = 'in', res = 800, pointsize = 8
# )
# netplot <- plot(
#   analysis.epithelial,
#   method = "diffusion",
#   data = fella.data,
#   nlimit = nlimit,
#   # LabelLengthAtPlot = 45,
#   NamesAsLabels = FALSE,
#   # font = 4,
#   plotLegend = T,
#   vertex.label.cex = .5,vertex.label.color="black"
#   # layout = TRUE,
#   # vertex.label.cex = vertex.label.cex
# )
# dev.off()
tab.all <- generateResultsTable(
  method = "diffusion",
  nlimit = nlimit,
  # LabelLengthAtPlot = 50,
  object = analysis.epithelial,
  data = fella.data
)
addWorksheet(
  wb, sheetName = l.a$subDir[i],
  gridLines = F, tabColour = '#008CCE'
)
writeDataTable(
  wb, sheet = i, tab.all, headerStyle = hs,
  tableStyle = 'TableStyleMedium27', withFilter = F
)
setColWidths(wb, sheet = i, cols = 1:ncol(tab.all), 'auto')
setRowHeights(wb, sheet = i, rows = 1, 24)
setRowHeights(wb, sheet = i, rows = 2:(nrow(tab.all) + 1), 16)
addStyle(
  wb, sheet = i, rows = 2:(nrow(tab.all) + 1),
  cols = 1:ncol(tab.all), gridExpand = T, stack = T,
  style = createStyle(valign = 'center')
)
}
saveWorkbook(
  wb, overwrite = T,
  paste0(
    l.a$output, '/results/',
    ifelse(is.na(pol), '', paste0(pol, '/')),
    'Network Analysis/',
    ifelse(is.na(pol), '', paste0(pol, '-')),
    'network.xlsx'
  )
)
}