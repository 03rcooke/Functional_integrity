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