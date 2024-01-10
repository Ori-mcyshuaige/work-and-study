library(vegan)

# Dissimilarity Indices for Community Ecologists
#data(varespec)
varespec<-read.table("./Biomarkers_Plasma_Urine_PERMANOVA_intensity.txt",header=TRUE,check=F,row.names=1)
varespec<-t(varespec)
# bray distance
vare.dist1 <- vegdist(varespec,diag = TRUE,upper = TRUE,method="bray")
# euclidean distance
vare.dist2 <- vegdist(varespec,diag = TRUE,upper = TRUE,method="euclidean")
vare.dist2=as.matrix(vare.dist2)
write.table(vare.dist2,"euclidean_distance.txt",sep='\t')



