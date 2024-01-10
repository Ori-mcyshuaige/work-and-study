library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(2408)
# load the Zhou et al snRNA-seq dataset
scRNA <- readRDS('final_scRNA_rpca.rds')
colnames(scRNA@meta.data)[4] <- 'tissue_type'
notUnknown <- subset(scRNA,subset = copykat.pred!='Unknown')


#####Set up Seurat object for WGCNA
notUnknown <- SetupForWGCNA(
  notUnknown,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "notUnknown" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
notUnknown <- MetacellsByGroups(
  seurat_obj = notUnknown,
  group.by = c('Cell.types','orig.ident'), # specify the columns in notUnknown@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 20, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  min_cells = 100,
  ident.group = 'Cell.types' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
notUnknown <- NormalizeMetacells(notUnknown)

#####Co-expression network analysis
notUnknown <- SetDatExpr(
  notUnknown,
  group_name = unique(notUnknown@misc$notUnknown$wgcna_metacell_obj$Cell.types), # the name of the group of interest in the group.by column
  group.by='Cell.types', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
notUnknown <- TestSoftPowers(
  notUnknown,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(notUnknown)

# assemble with patchwork
p <- wrap_plots(plot_list, ncol=2)
p
ggsave('05.hdwgcna/notUnknown_PlotSoftPowers.pdf',p,width = 7,height = 7)
#################Construct co-expression network
power_table <- GetPowerTable(notUnknown)
head(power_table)

# construct co-expression network:
notUnknown <- ConstructNetwork(
  notUnknown, soft_power=9,
  setDatExpr=FALSE
  # tom_name = 'all_cell' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(notUnknown, main='notUnknown hdWGCNA Dendrogram')

#########Module Eigengenes and Connectivity

# need to run ScaleData first or else harmony throws an error:
notUnknown <- ScaleData(notUnknown, features=VariableFeatures(notUnknown))

# compute all MEs in the full single-cell dataset
notUnknown <- ModuleEigengenes(
  notUnknown,
  group.by.vars="orig.ident"
)##

# harmonized module eigengenes:
hMEs <- GetMEs(notUnknown)

# module eigengenes:
MEs <- GetMEs(notUnknown, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
notUnknown <- ModuleConnectivity(
  notUnknown,
  group.by = 'Cell.types', group_name = unique(notUnknown@misc$notUnknown$wgcna_metacell_obj$Cell.types)
)

# rename the modules
# notUnknown <- ResetModuleNames(
#   notUnknown,
#   new_name = "hcc-cir"
# )

# plot genes ranked by kME for each module
p <- PlotKMEs(notUnknown, ncol=3)

p


# get the module assignment table:
modules <- GetModules(notUnknown)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(notUnknown, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
notUnknown <- ModuleExprScore(
  notUnknown,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
notUnknown <- ModuleExprScore(
  notUnknown,
  n_genes = 25,
  method='UCell'
)


######Basic Visualization


# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  notUnknown,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)


########Module Correlations
# plot module correlagram
ModuleCorrelogram(notUnknown)

# get hMEs from seurat object
MEs <- GetMEs(notUnknown, harmonized=TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
notUnknown@meta.data <- cbind(notUnknown@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(notUnknown, features=mods, group.by = 'Cell.types')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

# Plot INH-M4 hME using Seurat VlnPlot function
p <- VlnPlot(
  notUnknown,
  features = colnames(MEs),
  group.by = 'Cell.types',raster=FALSE,
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()

# plot output
p

# list of traits to correlate
cur_traits <- c('disulfidptosis_score', 'disulfidptosis_stage', 
                'tissue_type', 'malignant')

notUnknown$disulfidptosis_stage <- factor(notUnknown$disulfidptosis_stage,levels = c('Low','High'))
notUnknown$tissue_type <- factor(notUnknown$tissue_type,levels = c('cirrhotic','HCC'))
notUnknown@meta.data$malignant <- factor(notUnknown$copykat.pred,levels = c('non_Malignant','Malignant'))





notUnknown <- ModuleTraitCorrelation(
  notUnknown,
  traits = cur_traits,
  group.by='Cell.types'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(notUnknown)
names(mt_cor)
names(mt_cor$cor)
mt_cor$cor$all_cells

PlotModuleTraitCorrelation(
  notUnknown,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'black',
  high_color = 'red',
  mid_color = 'white',
  low_color = 'blue',
  plot_max = 0.5,
  combine=TRUE
)


########Hepatocytes

Hepatocytes <- subset(scRNA,subset = Cell.types=='Hepatocytes' & copykat.pred!='Unknown')


Hepatocytes <- SetupForWGCNA(
  Hepatocytes,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Hepatocytes" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
Hepatocytes <- MetacellsByGroups(
  seurat_obj = Hepatocytes,
  group.by = c('copykat.pred','orig.ident'), # specify the columns in Hepatocytes@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 7, # nearest-neighbors parameter
  max_shared = 3, # maximum number of shared cells between two metacells
  min_cells = 50,
  ident.group = 'copykat.pred' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Hepatocytes <- NormalizeMetacells(Hepatocytes)

#####Co-expression network analysis
Hepatocytes <- SetDatExpr(
  Hepatocytes,
  group_name = unique(Hepatocytes$copykat.pred), # the name of the group of interest in the group.by column
  group.by='copykat.pred', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
Hepatocytes <- TestSoftPowers(
  Hepatocytes,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(Hepatocytes)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

#################Construct co-expression network
power_table <- GetPowerTable(Hepatocytes)
head(power_table)

# construct co-expression network:
Hepatocytes <- ConstructNetwork(
  Hepatocytes, soft_power=8,overwrite_tom = TRUE,
  setDatExpr=FALSE
  # tom_name = 'all_cell' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(Hepatocytes, main='Hepatocytes hdWGCNA Dendrogram')

#########Module Eigengenes and Connectivity

# need to run ScaleData first or else harmony throws an error:
Hepatocytes <- ScaleData(Hepatocytes, features=VariableFeatures(Hepatocytes))

# compute all MEs in the full single-cell dataset
Hepatocytes <- ModuleEigengenes(
  Hepatocytes,
  group.by.vars="orig.ident"
)##

# harmonized module eigengenes:
hMEs <- GetMEs(Hepatocytes)

# module eigengenes:
MEs <- GetMEs(Hepatocytes, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Hepatocytes <- ModuleConnectivity(
  Hepatocytes,
  group.by = 'copykat.pred', group_name = unique(Hepatocytes$copykat.pred)
)

# rename the modules
# Hepatocytes <- ResetModuleNames(
#   Hepatocytes,
#   new_name = "hcc-cir"
# )

# plot genes ranked by kME for each module
p <- PlotKMEs(Hepatocytes, ncol=3)

p


# get the module assignment table:
modules <- GetModules(Hepatocytes)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(Hepatocytes, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
Hepatocytes <- ModuleExprScore(
  Hepatocytes,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
Hepatocytes <- ModuleExprScore(
  Hepatocytes,
  n_genes = 25,
  method='UCell'
)


######Basic Visualization


# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  Hepatocytes,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)


########Module Correlations
# plot module correlagram
ModuleCorrelogram(Hepatocytes)

# get hMEs from seurat object
MEs <- GetMEs(Hepatocytes, harmonized=TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
Hepatocytes@meta.data <- cbind(Hepatocytes@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(Hepatocytes, features=mods, group.by = 'copykat.pred')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

# Plot INH-M4 hME using Seurat VlnPlot function
p <- VlnPlot(
  Hepatocytes,
  features = colnames(MEs),
  group.by = 'copykat.pred',raster=FALSE,
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()

# plot output
p

# list of traits to correlate
cur_traits <- c('disulfidptosis_score', 'disulfidptosis_stage', 'tissue_type')
Hepatocytes$disulfidptosis_stage <- factor(Hepatocytes$disulfidptosis_stage,levels = c('Low','High'))
Hepatocytes$tissue_type <- factor(Hepatocytes$tissue_type,levels = c('cirrhotic','HCC'))


Hepatocytes <- ModuleTraitCorrelation(
  Hepatocytes,
  traits = cur_traits,
  group.by='copykat.pred'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(Hepatocytes)
names(mt_cor)
names(mt_cor$cor)
mt_cor$cor$all_cells

PlotModuleTraitCorrelation(
  Hepatocytes,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'black',
  high_color = 'red',
  mid_color = 'white',
  low_color = 'blue',
  plot_max = 0.5,
  combine=TRUE
)

