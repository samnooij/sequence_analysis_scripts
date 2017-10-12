#Open file

blast.output <- read.csv(file = "/data/sapo/experiments/20170411_pasc/results/blast/VP1_tblastn.bla", sep = "\t")
column.names <- c("query", "target", "identity", "length", "mismatches", "gap.openings", "q.start", "q.end", "t.start", "t.end", "evalue", "bit.score")
colnames(blast.output) <- column.names

blast.output <- within(blast.output, distance <- 100 - identity)

head(blast.output[,c("query", "target", "distance")])

blast.output$query.short <- sapply(strsplit(as.character(blast.output$query), "\\_"), "[", 1)
blast.output$query.short <- as.factor(blast.output$query.short)
blast.output$genogroup <- sapply(strsplit(as.character(blast.output$query.short), "[|]"), "[", 1)

head(blast.output[,c("query.short", "target", "distance")])
tail(blast.output[,c("query.short", "target", "distance")])

distance.df <- blast.output[,c("query.short", "target", "distance")]
identity.df <- blast.output[,c("query.short", "target", "identity")]

library(reshape2)

distance.mat <- acast(data = distance.df, formula = query.short ~ target, fun.aggregate = mean)
str(distance.mat)
identity.mat <- acast(data = identity.df, formula = query.short ~ target, fun.aggregate = mean)

##Apparently there are NAs in these matrices: replace them properly
distance.mat[8,8] <- 0
#there was only one here: distance to self was NA, is now 0
identity.mat[is.na(identity.mat)] <- 45
#BLAST didn't report anything lower than 50. Now those are replaced by 45.

bitscore.df <- blast.output[,c("query.short", "target", "bit.score")]
bitscore.mat <- acast(data = bitscore.df, formula = query.short ~ target, fun.aggregate = mean)

#### HEATMAP ####

library(pheatmap)

pheatmap(distance.mat)
pheatmap(identity.mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = c(10, 17, 30, 34, 36, 43, 45, 56, 62, 63, 64, 66, 68, 71),
  gaps_col = c(10, 17, 30, 34, 36, 43, 45, 56, 62, 63, 64, 66, 68, 71)
  )

pheatmap(identity.mat,
  clustering_method = "centroid")

pheatmap(identity.mat,
  clustering_method = "ward.D")

pheatmap(identity.mat,
  clustering_method = "ward.D2")

pheatmap(identity.mat,
  clustering_method = "median")


pheatmap(bitscore.mat, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  gaps_row = c(10, 17, 30, 34, 36, 43, 45, 56, 62, 63, 64, 66, 68, 71), 
  gaps_col = c(10, 17, 30, 34, 36, 43, 45, 56, 62, 63, 64, 66, 68, 71))

pheatmap(bitscore.mat,
  clustering_method = "ward.D")

pheatmap(bitscore.mat,
  clustering_method = "median")

#### FORCE-DIRECTED GRAPH ####
# https://stackoverflow.com/a/37832966

library(qgraph)

identity.mat <- acast(data = blast.output[,c("query.short", "target", "identity")], formula = query.short ~ target, fun.aggregate = mean)

qgraph(identity.mat, layout = "spring", vsize =3)
#nice visualisation, but probably works better when lines < [cutoff] identity are used (e.g. 50%)

identity.mat[ identity.mat<50 ] <- NA
qgraph(identity.mat, layout = "spring", vsize =3)
#that works reasonably well, but it is hard to read the labels...

#### CLUSTER DENDROGRAM ####
# http://www.statmethods.net/advstats/cluster.html

d <- dist(distance.mat, method = "euclidian") #another distance matrix?
fit <- hclust(d, method = "ward.D") #hierarchical clustering, like pheatmap
plot(fit)

groups <- cutree(fit, k = 5) #cut into k clusters
rect.hclust(fit, k=5, border="red") #and draw red borders around the clusters
#Above is just **an example**

#and with bootstrapped values:
library(pvclust)
fit <- pvclust(distance.mat, method.hclust = "ward.D",
  method.dist = "euclidian")
plot(fit)

#and put rectangles around highly supported groups:
pvrect(fit, alpha = .9)

### Model based clustering ###

library(mclust)

#Fit a model:
fit <- Mclust(data = distance.df)

#Plot the model results:
plot(fit)

summary(fit)
#So a VVV model (ellipsoidal with varying volume, shape and orientation) with 7 groups works best?
# I don't really get how to use this, even though it sounds interesting...

#### K-MEANS CLUSTERING ####

k.fit <- kmeans(x = distance.mat, centers = 7)

#Cluster plot: (2 principal components)
library(cluster)
clusplot(distance.mat, k.fit$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)
#Plots clusters in 2-dimensional space, drawing borders around k-means clusters

#Centroid plot: (2 discriminant functions)
library(fpc)
plotcluster(distance.mat, k.fit$cluster)
#plots centroids in 2-dimensional space, using hexadecimal numbers for centroids

## Validating cluster solutions, using fpc
k.fit.7 <- kmeans(x = distance.mat, centers = 7)
k.fit.15 <- kmeans(x = distance.mat, centers = 15)
cluster.stats(d = distance.mat, k.fit.7$cluster, k.fit.15$cluster)
cluster.stats(d = distance.mat, k.fit.15$cluster, k.fit.7$cluster)

# The above cluster.stats can be used to show number of clusters, members per cluster,
# distances between and within clusters and more.

#To calculate how much of the variance is explained by k-means, you have to divide
#the between sum of squares (betweenss) by the total sum of squares (totalss)

cluster.fun <- function(data, k) {
  cluster.kmeans <- kmeans(x = data, centers = k, iter.max = 10000, nstart = 30)
  cluster.list <- with(cluster.kmeans, list(cluster = cluster, var.explained = betweenss/totss))
  return(cluster.list)
}

#Find a good k-fit:
km.list <- lapply(X = 1:20, FUN = cluster.fun, data = distance.mat)
var.explained <- sapply(X = km.list, FUN = function(x) x$var.explained)

plot(1:20, var.explained, type = "b", xlab = "Number of clusters")
#Add lines to indicate "the best" point
abline(v = 8, h = var.explained[8], lty = 2, col = "grey")

#There seem to be bumps at 6, 8, 14 and 16
k6fit <- kmeans(x = distance.mat, centers = 6)
clusplot(distance.mat, k6fit$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0, main = "Kmeans 6 clusters")
k8fit <- kmeans(x = distance.mat, centers = 8)
clusplot(distance.mat, k8fit$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)
k14fit <- kmeans(x = distance.mat, centers = 14)
clusplot(distance.mat, k14fit$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)
k16fit <- kmeans(x = distance.mat, centers = 16)
clusplot(distance.mat, k16fit$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)

#show in a dendrogram how 8 clusters look:
d <- dist(distance.mat, method = "euclidian") #another distance matrix?
fit <- hclust(d, method = "ward.D") #hierarchical clustering, like pheatmap
plot(fit)

groups <- cutree(fit, k = 8) #cut into k clusters
rect.hclust(fit, k=8, border="red") #and draw red borders around the clusters

### Mutlidimensional scaling ###
# Metric MDS
d <- dist(distance.mat)
fit <- cmdscale(d,eig = TRUE, k=2)

x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS of VP1 tBLASTn", type = "n")
text(x, y, labels = row.names(distance.mat), cex = .7)
## ORIGINAL GENOGROUPS ##
text(x, y, labels = sapply(strsplit(row.names(distance.mat), "[|]"), "[", 1), cex = .7)
## till here ##

km.test <- cluster.fun(data = fit$points, k = 8)
plot(x = fit$points[,1], y = fit$points[,2], col = km.test$cluster)

plot.df <- as.data.frame(fit$points)
plot.df["cluster"] <- km.test$cluster

library(ggplot2)
col.clust.plot <- ggplot(data = plot.df, mapping = aes(x = V1, y = V2, colour = cluster)) +
  geom_point()
## ORIGINAL GENOGROUPS ##
library(RColorBrewer)
genogroup.palette <- c(brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 7, name = "Set3"))
genogroup <- sapply(strsplit(row.names(distance.mat), "[|]"), "[", 1)
col.clust.plot <- ggplot(data = plot.df, mapping = aes(x = V1, y = V2, colour = genogroup)) +
  geom_point(size = 3) +
  scale_colour_manual(values = genogroup.palette)
col.clust.plot
library(directlabels)
direct.label(col.clust.plot, method = "smart.grid")

plot(x = fit$points[,1], y = fit$points[,2], col = as.factor(genogroup), pch = 19)
text(x, y, labels = genogroup, cex = .5, pos = 3)


#Nonmetric MDS
library(MASS)
d <- dist(distance.mat)
fit.non <- isoMDS(d, k=2) #Doesn't work with zero distances
