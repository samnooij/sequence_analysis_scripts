#Cluster analysis, step by step

#1. import datafile:
blast.output <- read.csv(file = "/data/sapo/experiments/20170411_pasc/results/blast/VP1_tblastn.bla", sep = "\t")

#This is blast tabular output, label columns accordingly:
column.names <- c("query", "target", "identity", "length", "mismatches", "gap.openings", "q.start", "q.end", "t.start", "t.end", "evalue", "bit.score")
colnames(blast.output) <- column.names

#Make an additional column for short names (genogroup + accession ID - equal to "target")
blast.output$query.short <- sapply(strsplit(as.character(blast.output$query), "\\_"), "[", 1)
blast.output$query.short <- as.factor(blast.output$query.short)

#2. Make a column with distances (100 - % identity)
blast.output <- within(blast.output, distance <- 100 - identity)

#Make dataframes with only the required information
distance.df <- blast.output[,c("query.short", "target", "distance")]
identity.df <- blast.output[,c("query.short", "target", "identity")]

#And make matrices
library(reshape2)

distance.mat <- acast(data = distance.df, formula = query.short ~ target, fun.aggregate = mean)
identity.mat <- acast(data = identity.df, formula = query.short ~ target, fun.aggregate = mean)

#"Clean" the data by filling in NAs
distance.mat[8,8] <- 0
#there was only one here: distance to self was NA, is now 0
identity.mat[is.na(identity.mat)] <- 45
#BLAST didn't report anything lower than 50. Now those are replaced by 45.

#3. Quick look at the data with a **heatmap**

library(pheatmap)

pheatmap(distance.mat)
pheatmap(identity.mat, clustering_method = "centroid")

pheatmap(identity.mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = c(10, 17, 30, 34, 36, 43, 45, 56, 62, 63, 64, 66, 68, 71),
  gaps_col = c(10, 17, 30, 34, 36, 43, 45, 56, 62, 63, 64, 66, 68, 71)
)

#Calculate for each genogroup the within (intra) average distance
gi.avg <- mean(distance.mat[1:10,1:10])
gii.avg <- mean(distance.mat[11:17,11:17])
giii.avg <- mean(distance.mat[18:30,18:30])
giv.avg <- mean(distance.mat[31:34,31:34])
gv.avg <- mean(distance.mat[37:43,37:43])
gvi.avg <- mean(distance.mat[44:45,44:45])
gvii.avg <- mean(distance.mat[46:56,46:56])
gviii.avg <- mean(distance.mat[57:62,57:62])
gix.avg <- mean(distance.mat[35:36,35:36])
gx.avg <- distance.mat[63,63]
gxi.avg <- distance.mat[64,64]
gxii.avg <- mean(distance.mat[65:66,65:66])
gxiii.avg <- mean(distance.mat[67:68,67:68])
gxiv.avg <- mean(distance.mat[69:71,69:71])
gxv.avg <- mean(distance.mat[72:73,72:73])

#Calculate average distances between (inter) genogroups
gi.inter.gii <- (mean(distance.mat[1:10,11:17]) + mean(distance.mat[11:17,1:10])) / 2
gi.inter.giii <- (mean(distance.mat[1:10,18:30]) + mean(distance.mat[18:30,1:10])) / 2
gi.inter.giv <- (mean(distance.mat[1:10,31:34]) + mean(distance.mat[31:34,1:10])) / 2
gi.inter.gv <- (mean(distance.mat[1:10,37:43]) + mean(distance.mat[37:43,1:10])) / 2
gi.inter.gvi <- (mean(distance.mat[1:10,44:45]) + mean(distance.mat[44:45,1:10])) / 2
gi.inter.gvii <- (mean(distance.mat[1:10,46:56]) + mean(distance.mat[46:56,1:10])) / 2
gi.inter.gviii <- (mean(distance.mat[1:10,57:62]) + mean(distance.mat[57:62,1:10])) / 2
gi.inter.gix <- (mean(distance.mat[1:10,35:36]) + mean(distance.mat[35:36,1:10])) / 2
gi.inter.gx <- (mean(distance.mat[1:10,63]) + mean(distance.mat[63,1:10])) / 2
gi.inter.gxi <- (mean(distance.mat[1:10,64]) + mean(distance.mat[64,1:10])) / 2
gi.inter.gxii <- (mean(distance.mat[1:10,65:66]) + mean(distance.mat[65:66,1:10])) / 2
gi.inter.gxiii <- (mean(distance.mat[1:10,67:68]) + mean(distance.mat[67:68,1:10])) / 2
gi.inter.gxiv <- (mean(distance.mat[1:10,69:71]) + mean(distance.mat[69:71,1:10])) / 2
gi.inter.gxv <- (mean(distance.mat[1:10,72:73]) + mean(distance.mat[72:73,1:10])) / 2

gii.inter.giii <- (mean(distance.mat[11:17,18:30]) + mean(distance.mat[18:30,11:17])) / 2
gii.inter.giv <- (mean(distance.mat[11:17,31:34]) + mean(distance.mat[31:34,11:17])) / 2
gii.inter.gv <- (mean(distance.mat[11:17,37:43]) + mean(distance.mat[37:43,11:17])) / 2
gii.inter.gvi <- (mean(distance.mat[11:17,44:45]) + mean(distance.mat[44:45,11:17])) / 2
gii.inter.gvii <- (mean(distance.mat[11:17,46:56]) + mean(distance.mat[46:56,11:17])) / 2
gii.inter.gviii <- (mean(distance.mat[11:17,57:62]) + mean(distance.mat[57:62,11:17])) / 2
gii.inter.gix <- (mean(distance.mat[11:17,35:36]) + mean(distance.mat[35:36,11:17])) / 2
gii.inter.gx <- (mean(distance.mat[11:17,63]) + mean(distance.mat[63,11:17])) / 2
gii.inter.gxi <- (mean(distance.mat[11:17,64]) + mean(distance.mat[64,11:17])) / 2
gii.inter.gxii <- (mean(distance.mat[11:17,65:66]) + mean(distance.mat[65:66,11:17])) / 2
gii.inter.gxiii <- (mean(distance.mat[11:17,67:68]) + mean(distance.mat[67:68,11:17])) / 2
gii.inter.gxiv <- (mean(distance.mat[11:17,69:71]) + mean(distance.mat[69:71,11:17])) / 2
gii.inter.gxv <- (mean(distance.mat[11:17,72:73]) + mean(distance.mat[72:73,11:17])) / 2

giii.inter.giv <- (mean(distance.mat[18:30,31:34]) + mean(distance.mat[31:34,18:30])) / 2
giii.inter.gv <- (mean(distance.mat[18:30,37:43]) + mean(distance.mat[37:43,18:30])) / 2
giii.inter.gvi <- (mean(distance.mat[18:30,44:45]) + mean(distance.mat[44:45,18:30])) / 2
giii.inter.gvii <- (mean(distance.mat[18:30,46:56]) + mean(distance.mat[46:56,18:30])) / 2
giii.inter.gviii <- (mean(distance.mat[18:30,57:62]) + mean(distance.mat[57:62,18:30])) / 2
giii.inter.gix <- (mean(distance.mat[18:30,35:36]) + mean(distance.mat[35:36,18:30])) / 2
giii.inter.gx <- (mean(distance.mat[18:30,63]) + mean(distance.mat[63,18:30])) / 2
giii.inter.gxi <- (mean(distance.mat[18:30,64]) + mean(distance.mat[64,18:30])) / 2
giii.inter.gxii <- (mean(distance.mat[18:30,65:66]) + mean(distance.mat[65:66,18:30])) / 2
giii.inter.gxiii <- (mean(distance.mat[18:30,67:68]) + mean(distance.mat[67:68,18:30])) / 2
giii.inter.gxiv <- (mean(distance.mat[18:30,69:71]) + mean(distance.mat[69:71,18:30])) / 2
giii.inter.gxv <- (mean(distance.mat[18:30,72:73]) + mean(distance.mat[72:73,18:30])) / 2

giv.inter.gv <- (mean(distance.mat[31:34,37:43]) + mean(distance.mat[37:43,31:34])) / 2
giv.inter.gvi <- (mean(distance.mat[31:34,44:45]) + mean(distance.mat[44:45,31:34])) / 2
giv.inter.gvii <- (mean(distance.mat[31:34,46:56]) + mean(distance.mat[46:56,31:34])) / 2
giv.inter.gviii <- (mean(distance.mat[31:34,57:62]) + mean(distance.mat[57:62,31:34])) / 2
giv.inter.gix <- (mean(distance.mat[31:34,35:36]) + mean(distance.mat[35:36,31:34])) / 2
giv.inter.gx <- (mean(distance.mat[31:34,63]) + mean(distance.mat[63,31:34])) / 2
giv.inter.gxi <- (mean(distance.mat[31:34,64]) + mean(distance.mat[64,31:34])) / 2
giv.inter.gxii <- (mean(distance.mat[31:34,65:66]) + mean(distance.mat[65:66,31:34])) / 2
giv.inter.gxiii <- (mean(distance.mat[31:34,67:68]) + mean(distance.mat[67:68,31:34])) / 2
giv.inter.gxiv <- (mean(distance.mat[31:34,69:71]) + mean(distance.mat[69:71,31:34])) / 2
giv.inter.gxv <- (mean(distance.mat[31:34,72:73]) + mean(distance.mat[72:73,31:34])) / 2

gv.inter.gvi <- (mean(distance.mat[37:43,44:45]) + mean(distance.mat[44:45,37:43])) / 2
gv.inter.gvii <- (mean(distance.mat[37:43,46:56]) + mean(distance.mat[46:56,37:43])) / 2
gv.inter.gviii <- (mean(distance.mat[37:43,57:62]) + mean(distance.mat[57:62,37:43])) / 2
gv.inter.gix <- (mean(distance.mat[37:43,35:36]) + mean(distance.mat[35:36,37:43])) / 2
gv.inter.gx <- (mean(distance.mat[37:43,63]) + mean(distance.mat[63,37:43])) / 2
gv.inter.gxi <- (mean(distance.mat[37:43,64]) + mean(distance.mat[64,37:43])) / 2
gv.inter.gxii <- (mean(distance.mat[37:43,65:66]) + mean(distance.mat[65:66,37:43])) / 2
gv.inter.gxiii <- (mean(distance.mat[37:43,67:68]) + mean(distance.mat[67:68,37:43])) / 2
gv.inter.gxiv <- (mean(distance.mat[37:43,69:71]) + mean(distance.mat[69:71,37:43])) / 2
gv.inter.gxv <- (mean(distance.mat[37:43,72:73]) + mean(distance.mat[72:73,37:43])) / 2

gvi.inter.gvii <- (mean(distance.mat[44:45,46:56]) + mean(distance.mat[46:56,44:45])) / 2
gvi.inter.gviii <- (mean(distance.mat[44:45,57:62]) + mean(distance.mat[57:62,44:45])) / 2
gvi.inter.gix <- (mean(distance.mat[44:45,35:36]) + mean(distance.mat[35:36,44:45])) / 2
gvi.inter.gx <- (mean(distance.mat[44:45,63]) + mean(distance.mat[63,44:45])) / 2
gvi.inter.gxi <- (mean(distance.mat[44:45,64]) + mean(distance.mat[64,44:45])) / 2
gvi.inter.gxii <- (mean(distance.mat[44:45,65:66]) + mean(distance.mat[65:66,44:45])) / 2
gvi.inter.gxiii <- (mean(distance.mat[44:45,67:68]) + mean(distance.mat[67:68,44:45])) / 2
gvi.inter.gxiv <- (mean(distance.mat[44:45,69:71]) + mean(distance.mat[69:71,44:45])) / 2
gvi.inter.gxv <- (mean(distance.mat[44:45,72:73]) + mean(distance.mat[72:73,44:45])) / 2

gvii.inter.gviii <- (mean(distance.mat[46:56,57:62]) + mean(distance.mat[57:62,46:56])) / 2
gvii.inter.gix <- (mean(distance.mat[46:56,35:36]) + mean(distance.mat[35:36,46:56])) / 2
gvii.inter.gx <- (mean(distance.mat[46:56,63]) + mean(distance.mat[63,46:56])) / 2
gvii.inter.gxi <- (mean(distance.mat[46:56,64]) + mean(distance.mat[64,46:56])) / 2
gvii.inter.gxii <- (mean(distance.mat[46:56,65:66]) + mean(distance.mat[65:66,46:56])) / 2
gvii.inter.gxiii <- (mean(distance.mat[46:56,67:68]) + mean(distance.mat[67:68,46:56])) / 2
gvii.inter.gxiv <- (mean(distance.mat[46:56,69:71]) + mean(distance.mat[69:71,46:56])) / 2
gvii.inter.gxv <- (mean(distance.mat[46:56,72:73]) + mean(distance.mat[72:73,46:56])) / 2

gviii.inter.gix <- (mean(distance.mat[57:62,35:36]) + mean(distance.mat[35:36,57:62])) / 2
gviii.inter.gx <- (mean(distance.mat[57:62,63]) + mean(distance.mat[63,57:62])) / 2
gviii.inter.gxi <- (mean(distance.mat[57:62,64]) + mean(distance.mat[64,57:62])) / 2
gviii.inter.gxii <- (mean(distance.mat[57:62,65:66]) + mean(distance.mat[65:66,57:62])) / 2
gviii.inter.gxiii <- (mean(distance.mat[57:62,67:68]) + mean(distance.mat[67:68,57:62])) / 2
gviii.inter.gxiv <- (mean(distance.mat[57:62,69:71]) + mean(distance.mat[69:71,57:62])) / 2
gviii.inter.gxv <- (mean(distance.mat[57:62,72:73]) + mean(distance.mat[72:73,57:62])) / 2

gix.inter.gx <- (mean(distance.mat[35:36,63]) + mean(distance.mat[63,35:36])) / 2
gix.inter.gxi <- (mean(distance.mat[35:36,64]) + mean(distance.mat[64,35:36])) / 2
gix.inter.gxii <- (mean(distance.mat[35:36,65:66]) + mean(distance.mat[65:66,35:36])) / 2
gix.inter.gxiii <- (mean(distance.mat[35:36,67:68]) + mean(distance.mat[67:68,35:36])) / 2
gix.inter.gxiv <- (mean(distance.mat[35:36,69:71]) + mean(distance.mat[69:71,35:36])) / 2
gix.inter.gxv <- (mean(distance.mat[35:36,72:73]) + mean(distance.mat[72:73,35:36])) / 2

gx.inter.gxi <- (mean(distance.mat[63,64]) + mean(distance.mat[64,63])) / 2
gx.inter.gxii <- (mean(distance.mat[63,65:66]) + mean(distance.mat[65:66,63])) / 2
gx.inter.gxiii <- (mean(distance.mat[63,67:68]) + mean(distance.mat[67:68,63])) / 2
gx.inter.gxiv <- (mean(distance.mat[63,69:71]) + mean(distance.mat[69:71,63])) / 2
gx.inter.gxv <- (mean(distance.mat[63,72:73]) + mean(distance.mat[72:73,63])) / 2

gxi.inter.gxii <- (mean(distance.mat[64,65:66]) + mean(distance.mat[65:66,64])) / 2
gxi.inter.gxiii <- (mean(distance.mat[64,67:68]) + mean(distance.mat[67:68,64])) / 2
gxi.inter.gxiv <- (mean(distance.mat[64,69:71]) + mean(distance.mat[69:71,64])) / 2
gxi.inter.gxv <- (mean(distance.mat[64,72:73]) + mean(distance.mat[72:73,64])) / 2

gxii.inter.gxiii <- (mean(distance.mat[65:66,67:68]) + mean(distance.mat[67:68,65:66])) / 2
gxii.inter.gxiv <- (mean(distance.mat[65:66,69:71]) + mean(distance.mat[69:71,65:66])) / 2
gxii.inter.gxv <- (mean(distance.mat[65:66,72:73]) + mean(distance.mat[72:73,65:66])) / 2

gxiii.inter.gxiv <- (mean(distance.mat[67:68,69:71]) + mean(distance.mat[69:71,67:68])) / 2
gxiii.inter.gxv <- (mean(distance.mat[67:68,72:73]) + mean(distance.mat[72:73,67:68])) / 2

gxiv.inter.gxv <- (mean(distance.mat[69:71,72:73]) + mean(distance.mat[72:73,69:71])) / 2


### Put all the inter- and intra-genogroup distances in a matrix: ###
inter.intra.matrix <- matrix(nrow = 15, ncol = 15)

inter.intra.matrix[1,] <- c(gi.avg, gi.inter.gii, gi.inter.giii, 
  gi.inter.giv, gi.inter.gv, gi.inter.gvi, gi.inter.gvii, gi.inter.gviii, 
  gi.inter.gix, gi.inter.gx, gi.inter.gxi, gi.inter.gxii, gi.inter.gxiii,
  gi.inter.gxiv, gi.inter.gxv)
inter.intra.matrix[2,2:15] <- c(gii.avg, gii.inter.giii, gii.inter.giv, 
  gii.inter.gv, gii.inter.gvi, gii.inter.gvii, gii.inter.gviii, gii.inter.gix, 
  gii.inter.gx, gii.inter.gxi, gii.inter.gxii, gii.inter.gxiii, gii.inter.gxiv, 
  gii.inter.gxv)
inter.intra.matrix[3,3:15] <- c(giii.avg, giii.inter.giv, giii.inter.gv,
  giii.inter.gvi, giii.inter.gvii, giii.inter.gviii, giii.inter.gix,
  giii.inter.gx, giii.inter.gxi, giii.inter.gxii, giii.inter.gxiii,
  giii.inter.gxiv, giii.inter.gxv)
inter.intra.matrix[4,4:15] <- c(giv.avg, giv.inter.gv, giv.inter.gvi,
  giv.inter.gvii, giv.inter.gviii, giv.inter.gix, giv.inter.gx,
  giv.inter.gxi, giv.inter.gxii, giv.inter.gxiii, giv.inter.gxiv, giv.inter.gxv)
inter.intra.matrix[5,5:15] <- c(gv.avg, gv.inter.gvi, gv.inter.gvii, 
  gv.inter.gviii, gv.inter.gix, gv.inter.gx, gv.inter.gxi, gv.inter.gxii,
  gv.inter.gxiii, gv.inter.gxiv, gv.inter.gxv)
inter.intra.matrix[6,6:15] <- c(gvi.avg, gvi.inter.gvii, gvi.inter.gviii,
  gvi.inter.gix, gvi.inter.gx, gvi.inter.gxi, gvi.inter.gxii, gvi.inter.gxiii,
  gvi.inter.gxiv, gvi.inter.gxv)
inter.intra.matrix[7,7:15] <- c(gvii.avg, gvii.inter.gviii, gvii.inter.gix, 
  gvii.inter.gx, gvii.inter.gxi, gvii.inter.gxii, gvii.inter.gxiii,
  gvii.inter.gxiv, gvii.inter.gxv)
inter.intra.matrix[8,8:15] <- c(gviii.avg, gviii.inter.gix, gviii.inter.gx,
  gviii.inter.gxi, gviii.inter.gxii, gviii.inter.gxiii, gviii.inter.gxiv,
  gviii.inter.gxv)
inter.intra.matrix[9,9:15] <- c(gix.avg, gix.inter.gx, gix.inter.gxi,
  gix.inter.gxii, gix.inter.gxiii, gix.inter.gxiv, gix.inter.gxv)
inter.intra.matrix[10,10:15] <- c(gx.avg, gx.inter.gxi, gx.inter.gxii,
  gx.inter.gxiii, gx.inter.gxiv, gx.inter.gxv)
inter.intra.matrix[11,11:15] <- c(gxi.avg, gxi.inter.gxii, gxi.inter.gxiii,
  gxi.inter.gxiv, gxi.inter.gxv)
inter.intra.matrix[12,12:15] <- c(gxii.avg, gxii.inter.gxiii, gxii.inter.gxiv,
  gxii.inter.gxv)
inter.intra.matrix[13,13:15] <- c(gxiii.avg, gxiii.inter.gxiv, gxiii.inter.gxv)
inter.intra.matrix[14,14:15] <- c(gxiv.avg, gxiv.inter.gxv)
inter.intra.matrix[15,15] <- c(gxv.avg)

colnames(inter.intra.matrix) <- c("GI", "GII", "GIII", "GIV", "GV", "GVI", "GVII", "GVIII", "GIX", "GX", "GXI", "GXII", "GXIII", "GXIV", "GXV")
rownames(inter.intra.matrix) <- c("GI", "GII", "GIII", "GIV", "GV", "GVI", "GVII", "GVIII", "GIX", "GX", "GXI", "GXII", "GXIII", "GXIV", "GXV")

### Plot the inter/intra matrix as heatmap: ###

pheatmap(inter.intra.matrix, 
  cluster_cols = FALSE, 
  cluster_rows = FALSE, 
  main = "Inter- and intra-group average distances", 
  display_numbers = TRUE)

#Calculate for each genogroup the average distance to all others
gi.inter <- (mean(distance.mat[1:10,11:73]) +
    mean(distance.mat[11:73,1:10])) / 2
gii.inter <- (mean(distance.mat[1:10,11:17]) +
    mean(distance.mat[11:17,1:10]) + 
    mean(distance.mat[19:73,11:18]) + 
    mean(distance.mat[11:18,19:73])) / 4
giii.inter <- (mean(distance.mat[1:17,18:30]) +
    mean(distance.mat[18:30,1:17]) +
    mean(distance.mat[31:73,18:30]) +
    mean(distance.mat[18:30,31:73])) / 4
giv.inter <- (mean(distance.mat[1:30,31:34]) +
    mean(distance.mat[31:34,1:30]) +
    mean(distance.mat[35:73,31:34]) +
    mean(distance.mat[31:34,35:73])) / 4
gv.inter <- (mean(distance.mat[1:36,37:43]) +
    mean(distance.mat[37:43,1:36]) +
    mean(distance.mat[44:73,37:43]) +
    mean(distance.mat[37:43,44:73])) / 4
gvi.inter <- (mean(distance.mat[1:43,44:45]) +
    mean(distance.mat[44:45,1:43]) +
    mean(distance.mat[46:73,44:45]) +
    mean(distance.mat[44:45,46:73])) / 4
gvii.inter <- (gi.inter.gvii + gii.inter.gvii + giii.inter.gvii + giv.inter.gvii + gv.inter.gvii + gvi.inter.gvii + gvii.inter.gviii + gvii.inter.gix + gvii.inter.gx + gvii.inter.gxi + gvii.inter.gxii + gvii.inter.gxiii + gvii.inter.gxiv + gvii.inter.gxv) / 14
gviii.inter <- (gi.inter.gviii + gii.inter.gviii + giii.inter.gviii + giv.inter.gviii + gv.inter.gviii + gvi.inter.gviii + gvii.inter.gviii + gviii.inter.gix + gviii.inter.gx + gviii.inter.gxi + gviii.inter.gxii + gviii.inter.gxiii + gviii.inter.gxiv + gviii.inter.gxv) / 14
gix.inter <- (gi.inter.gix + gii.inter.gix + giii.inter.gix + giv.inter.gix + gv.inter.gix + gvi.inter.gix + gvii.inter.gix + gviii.inter.gix + gix.inter.gx + gix.inter.gxi + gix.inter.gxii + gix.inter.gxiii + gix.inter.gxiv + gix.inter.gxv) / 14 
gx.inter <- (gi.inter.gx + gii.inter.gx + giii.inter.gx + giv.inter.gx + gv.inter.gx + gvi.inter.gx + gvii.inter.gx + gviii.inter.gx + gix.inter.gx + gx.inter.gxi + gx.inter.gxii + gx.inter.gxiii + gx.inter.gxiv + gx.inter.gxv) / 14
gxi.inter <- (gi.inter.gxi + gii.inter.gxi + giii.inter.gxi + giv.inter.gxi + gv.inter.gxi + gvi.inter.gxi + gvii.inter.gxi + gviii.inter.gxi + gix.inter.gxi + gx.inter.gxi + gxi.inter.gxii + gxi.inter.gxiii + gxi.inter.gxiv + gxi.inter.gxv) / 14
gxii.inter <- (gi.inter.gxii + gii.inter.gxii + giii.inter.gxii + giv.inter.gxii + gv.inter.gxii + gvi.inter.gxii + gvii.inter.gxii + gviii.inter.gxii + gix.inter.gxii + gx.inter.gxii + gxi.inter.gxii + gxii.inter.gxiii + gxii.inter.gxiv + gxii.inter.gxv) / 14
gxiii.inter <- (gi.inter.gxiii + gii.inter.gxiii + giii.inter.gxiii + giv.inter.gxiii + gv.inter.gxiii + gvi.inter.gxiii + gvii.inter.gxiii + gviii.inter.gxiii + gix.inter.gxiii + gx.inter.gxiii + gxi.inter.gxiii + gxii.inter.gxiii + gxiii.inter.gxiv + gxiii.inter.gxv) / 14
gxiv.inter <- (gi.inter.gxiv + gii.inter.gxiv + giii.inter.gxiv + giv.inter.gxiv + gv.inter.gxiv + gvi.inter.gxiv + gvii.inter.gxiv + gviii.inter.gxiv + gix.inter.gxiv + gx.inter.gxiv + gxi.inter.gxiv + gxii.inter.gxiv + gxiii.inter.gxiv + gxiv.inter.gxv) / 14
gxv.inter <- (gi.inter.gxv + gii.inter.gxv + giii.inter.gxv + giv.inter.gxv + gv.inter.gxv + gvi.inter.gxv + gvii.inter.gxv + gviii.inter.gxv + gix.inter.gxv + gx.inter.gxv + gxi.inter.gxv + gxii.inter.gxv + gxiii.inter.gxv + gxiv.inter.gxv) / 14

overall.mean.inter.distance <- (gi.inter + gii.inter + giii.inter +
    giv.inter + gv.inter + gvi.inter +
    gvii.inter + gviii.inter + gix.inter +
    gx.inter + gxi.inter + gxii.inter +
    gxiii.inter + gxiv.inter + gxv.inter) / 15
overall.mean.distance <- mean(distance.mat)

#4. Plot as **histograms**
library(ggplot2)
hist1 <- ggplot(data = identity.df, mapping = aes(x=identity)) + geom_histogram(breaks = seq(from = ceiling(min(identity.df$identity)), to = 100, by = 1), fill = "cyan", colour = "cyan4") + geom_rug() + labs(title = "All vs. all tBLASTn identity scores for Sapovirus VP1 genes", x = "Identity (%)", y = "Frequency")
hist1.line <- hist1 + geom_density(fill = "gold", alpha = .3, aes(y=..count..))
hist1.line

### CLUSTER ANALYSIS LIKE WITH UCLUST, WITH IDENTITY CUTOFFS ###
avg.clusters <- hclust(dist(distance.mat), method = "average") #average=UPGMA
plot(avg.clusters)
# distance cutoff:
avg.cut <- cutree(avg.clusters, h = 67)
number.of.clusters <- max(avg.cut)
# Either one of these visualisations:
plot(avg.clusters, labels = as.character(avg.cut))
rect.hclust(avg.clusters, h = 67)

#Save the results into a new dataframe:
cluster.df <- as.data.frame(avg.cut)
cluster.df$name <- rownames(cluster.df)
distance.df.clus <- merge.data.frame(x = distance.df, y = cluster.df, by.x = "query.short", by.y = "name")
colnames(distance.df.clus)[colnames(distance.df.clus)=="avg.cut"] <- "clusternumber"

cluster.df2 <- aggregate(as.character(query.short) ~ clusternumber, distance.df.clus, unique)
colnames(cluster.df2) <- c("cluster.id", "member.names")
# Count number of members per cluster:
cluster.df2$no.members <- sapply(cluster.df2[,2], FUN = length)
# the "member" column needs to be converted from list to character (comma separated)
#  to keep the number of dimensions in the dataframe right:
cluster.df2$member.names <- sapply(cluster.df2$member.names, paste, collapse = ", ")

#Export this dataframe:
setwd("/data/sapo/experiments/20170411_pasc/")
write.table(x = cluster.df2[c("cluster.id", "no.members", "member.names")], file = "tmp/distCutoff-Clusters_and_members.csv", sep = "\t", row.names = FALSE)

#Now do the 'cluster analysis' and plot (cutoff values - number of clusters)
cutoff.max = 150

#Record the numbers in a dataframe (via matrix)
cluster.analysis.mat <- matrix(ncol = 2, nrow = cutoff.max)

for (cutoff in cutoff.max:0) {
  tree.cut <- cutree(democlust, h = cutoff)
  number.of.clusters <- max(tree.cut)
  print(paste(cutoff, number.of.clusters, sep = " "))
  cluster.analysis.mat[cutoff,1] <- cutoff
  cluster.analysis.mat[cutoff,2] <- number.of.clusters
}

#Convert matrix to dataframe (can be exported)
cluster.analysis.df <- data.frame(cluster.analysis.mat)
colnames(cluster.analysis.df) <- c("cutoff", "number.of.clusters")

#Find the plateaus in number of clusters (5 consecutive numbers)
find.consecutives <- function(x, k = 5) {
  #Find the length of vector x
  n <- length(x)
  #Create an empty object for end points
  end.points <- NULL
  #Loop over x, from start+5 to k
  for (i in k:n) {
    #if 5 or more elements before i are equal to i, save the endpoint
    if (all(x[(i-k):i] == x[i])) end.points <- c(end.points, i)
  }
  new.points <- NULL
  for (i in 1:(length(end.points) -1)) {
    if (!end.points[i+1] == (end.points[i] + 1)) new.points <- c(new.points, end.points[i])
  }
  #Return endpoints
  return(new.points)
}

cutoff.intersects <- find.consecutives(cluster.analysis.df$number.of.clusters)

#Plot the cluster analysis results
plot(cluster.analysis.df, type = "b", xlab = "Cutoff value", ylab = "Number of clusters", main = "Cluster analysis for Sapovirus VP1 tBLASTn results")

#And record the results in another dataframe (via matrix)
intersect.mat <- matrix(ncol = 2, nrow = length(cutoff.intersects))

#Difficult way to plot the ablines and fill a dataframe with all values at the same time:
for (value in cutoff.intersects) {
  abline(v = value, h = cluster.analysis.df[value,2], lty = 2, col = "grey")
  intersect.mat[which(cutoff.intersects == value, arr.ind = TRUE), 1] = value
  intersect.mat[which(cutoff.intersects == value, arr.ind = TRUE), 2] = cluster.analysis.df$number.of.clusters[value]
}

#Matrix to dataframe (can be exported)
intersect.df <- data.frame(intersect.mat)
colnames(intersect.df) <- c("cutoff", "number.of.clusters")

intersect.df

"""
democlust2 <- hclust(dist(distance.mat), method = "mcquitty") #mcquitty=WPGMA
plot(democlust2)
rect.hclust(democlust2, h = 70)

demo3 <-  hclust(dist(distance.mat), method = "centroid") #centroid=UPGMC
plot(demo3)
rect.hclust(demo3, h = 55)
demo4 <- hclust(dist(distance.mat), method = "single")
plot(demo4)
rect.hclust(demo4, h = 64)
#64 does not separate GXIII-GV-GII
rect.hclust(demo4, h = 60)
#60 separates them better, but also cuts up GVI in 3 clusters
demo5 <- hclust(dist(distance.mat), method = "ward.D2")
plot(demo5)
rect.hclust(demo5, h = 110)
#110 does not separate GII from GXIII and GIX from GX,
# however, it cuts GIII and GVII into 2 groups each
demo6 <- hclust(dist(distance.mat), method = "complete")
plot(demo6)
rect.hclust(demo6, h = 70)
#70 does not separate GII-GV, while it cuts GI, GIII and GVII into 2 groups
demo7 <- hclust(dist(distance.mat), method = "median") #median=WPGMC
plot(demo7)
rect.hclust(demo7, h = 48)
#48 is too strict: GI, GIX & GXIII into 2, GVIII into 3
# but 49 is too lenient...: GIX+GVII, GVIII+GV+GII+GIII
"""

# x. k-means clustering
#Find the best number of clusters:
cluster.fun <- function(data, k) {
  cluster.kmeans <- kmeans(x = data, centers = k, iter.max = 10000, nstart = 30)
  cluster.list <- with(cluster.kmeans, list(cluster = cluster, var.explained = betweenss/totss))
  return(cluster.list)
}

#Find a **good k-fit**:
km.list <- lapply(X = 1:20, FUN = cluster.fun, data = distance.mat)
var.explained <- sapply(X = km.list, FUN = function(x) x$var.explained)

plot(1:20, var.explained, type = "b", xlab = "Number of clusters", ylab = "Variation explained (ratio ss)", main = "k-means analysis of Sapovirus VP1 tBLASTn results")

#Find an optimum by searching for the point where the variation grows less than 5%
k.means.df <- data.frame(1:20, var.explained)
colnames(k.means.df) <- c("number.of.clusters", "var.explained")

find.level.off <- function(x, l = 5) {
  #Find the length of vector x
  n <- length(x)
  #Loop over x
  for (i in 2:n) {
    #if the next element is less than 5% greater, the line has levelled off.
    if (x[i] <= (x[i-1] * 1.05)) return(x[i])
  }
}

level.off <- find.level.off(k.means.df$var.explained)
n.clusters <- which(k.means.df$var.explained == level.off, arr.ind = TRUE)
#Add lines to indicate "the best" point
abline(v = n.clusters, h = level.off, lty = 2, col = "grey")

#Draw a **Cluster Dendrogram** illustrating the fitting k-means
d <- dist(distance.mat, method = "euclidian") #another distance matrix?
fit <- hclust(d, method = "ward.D2") #hierarchical clustering, like pheatmap
plot(fit)

groups <- cutree(fit, k = n.clusters) #cut into k clusters
rect.hclust(fit, k=n.clusters, border="red") #and draw red borders around the clusters

#Write the results into a dataframe
cluster.k.df <- as.data.frame(groups)
cluster.k.df$name <- rownames(cluster.k.df)
colnames(cluster.k.df) <- c("cluster.id", "member.names")
cluster.k.df <- aggregate(member.names ~ cluster.id, cluster.k.df, unique)
# Count number of members per cluster:
cluster.k.df$no.members <- sapply(cluster.k.df[,2], FUN = length)

#And export the dataframe to csv
cluster.k.df$member.names <- sapply(cluster.k.df$member.names, paste, collapse = ", ")
write.table(x = cluster.k.df[c("cluster.id", "no.members", "member.names")], file = "tmp/k-Clusters_and_members.csv", sep = "\t", row.names = FALSE)

### Display/export **table** telling **clustering statistics**
library(fpc)
k.fit.9 <- kmeans(x = distance.mat, centers = 9)
clust.st <- cluster.stats(d = distance.mat, clustering = k.fit.9$cluster)

#Cluster sizes:
clust.st$cluster.size

# Average distances:
clust.st$average.distance

# "Separation":
clust.st$separation

#Average between clusters:
clust.st$average.between

#Average within clusters:
clust.st$average.within

#Dunn index (1 & 2) - note that these are not very high = not good, but with k = 15 it's much lower
c(clust.st$dunn, clust.st$dunn2)

### As an extra: **bootstrapped cluster dendrogram**:
library(pvclust)
fit <- pvclust(distance.mat, method.hclust = "ward.D2",
  method.dist = "euclidian")
plot(fit)

#and put rectangles around highly supported groups:
pvrect(fit, alpha = .90)
#The higher the cutoff, the more groups you get (for 70-95%)