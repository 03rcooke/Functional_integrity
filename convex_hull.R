### Plotting function to plot convex hulls
### Filename: Plot_ConvexHull.R
### Notes:
############################################################################

# INPUTS:
# xcoords: x-coordinates of point data
# ycoords: y-coordinates of point data
# lcolor: line color

# OUTPUTS:
# convex hull around data points in a particular color (specified by lcolor)

# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  
# END OF FUNCTION

library(ggplot2)

df     <- data.frame(iris)                   # iris dataset
pca    <- prcomp(df[,1:4], retx=T, scale.=T) # scaled pca [exclude species col]
scores <- pca$x[,1:3]                        # scores for first three PC's

# k-means clustering [assume 3 clusters]
km     <- kmeans(scores, centers=3, nstart=5)
ggdata <- data.frame(scores, Cluster=km$cluster, Species=df$Species)

# stat_ellipse is not part of the base ggplot package
source("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R") 

ggplot(ggdata) +
  geom_point(aes(x=PC1, y=PC2, color=factor(Cluster)), size=5, shape=20) +
  stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster"))














##################### Lefcheck ###################

#Visualize species' differences in multivariate trait space
#Perform k-means clustering with no a priori specification for k
traits_k <- pamk(gd, krange = 2:10)
#Perform multidimensional scaling on functional dendrogram
traits_nmds <- metaMDS(gd, k = c_, trymax = 500)
#Plot in two dimensions
#par(mar=c(4,4,1,1))
ordiplot(traits_nmds, type = "n")
# Assign colors to different groups
groups <- levels(factor(traits_k$pamobject$clustering))
points.symbols=15:16
points.colors=c("firebrick3","cornflowerblue")
for(i in seq_along(groups)) {
  points(traits_nmds$points[traits_k$pamobject$clustering==groups[i],],
         pch=points.symbols[i],col=points.colors[i],cex=1.4) }
ordispider(traits_nmds,factor(traits_k$pamobject$clustering),label=F)
ordihull(traits_nmds,factor(traits_k$pamobject$clustering),lty="dotted")
orditorp(traits_nmds,dis="sites",pcex=0,air=0.5,col="grey10",cex=0.8)


