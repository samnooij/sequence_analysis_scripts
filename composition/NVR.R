# Custom implementation of NVR (Natural Vector)
# Yu et al., 2013. PLoS ONE, Real Time Classification of Viruses in 12 Dimensions
# DOI: 10.1371/journal.pone.0064328

# Author: Sam Nooij
# Date: 10 Oct 2017
# Email: sam.nooij@rivm.nl

library('seqinr')
library('data.table')

nucleotides <- c('a', 'c', 'g', 't')

sequences <- read.fasta(file = "sequences.fasta")

nvr.df <- data.frame(
  name = getName(sequences))

frequency.df <- as.data.frame(do.call(rbind, lapply(X = sequences, FUN = table)))
frequency.df$name <- rownames(frequency.df)
setnames(frequency.df, old = c('a', 'c', 'g', 't'),
  new = c("a.frequency", "c.frequency", "g.frequency", "t.frequency"))

nvr.df <- merge(x = nvr.df, y = frequency.df)

mean.nucl.positions <- function(sequence){
  mean.a.position <- sum(which(sequence %in% 'a') / table(sequence)['a'])
  mean.c.position <- sum(which(sequence %in% 'c') / table(sequence)['c'])
  mean.g.position <- sum(which(sequence %in% 'g') / table(sequence)['g'])
  mean.t.position <- sum(which(sequence %in% 't') / table(sequence)['t'])
  return(c(mean.a.position, mean.c.position, mean.g.position, mean.t.position))
}

mean.pos.df <- as.data.frame(do.call(rbind, lapply(X = sequences, FUN = mean.nucl.positions)))
mean.pos.df$name <- rownames(mean.pos.df)
setnames(mean.pos.df, old = c('V1', 'V2', 'V3', 'V4'),
  new = c("mean.a.position", "mean.c.position", "mean.g.position", "mean.t.position"))

nvr.df <- merge(nvr.df, mean.pos.df)

central.moment <- function(sequence){
  a.positions <- which(sequence %in% 'a')
  c.positions <- which(sequence %in% 'c')
  g.positions <- which(sequence %in% 'g')
  t.positions <- which(sequence %in% 't')
  mean.a.position <- sum(a.positions / length(a.positions))
  mean.c.position <- sum(c.positions / length(c.positions))
  mean.g.position <- sum(g.positions / length(g.positions))
  mean.t.position <- sum(t.positions / length(t.positions))
  central.moment.a <- sum((a.positions - mean.a.position) ^ 2 / (length(a.positions) * length(sequence)))
  central.moment.c <- sum((c.positions - mean.c.position) ^ 2 / (length(c.positions) * length(sequence)))
  central.moment.g <- sum((g.positions - mean.g.position) ^ 2 / (length(g.positions) * length(sequence)))
  central.moment.t <- sum((t.positions - mean.t.position) ^ 2 / (length(t.positions) * length(sequence)))
  return(c(central.moment.a, central.moment.c, central.moment.g, central.moment.t))
}

central.moment.df <- as.data.frame(do.call(rbind, lapply(X = sequences, FUN = central.moment)))
central.moment.df$name <- rownames(central.moment.df)
setnames(central.moment.df, old = c('V1', 'V2', 'V3', 'V4'),
  new = c("central.moment.a", "central.moment.c", "central.moment.g", "central.moment.t"))

nvr.df <- merge(nvr.df, central.moment.df)

write.table(x = nvr.df, file = "sequences.nvr", sep = ",", row.names = FALSE, col.names = TRUE)

distance.matrix <- dist(x = nvr.df, method = "euclidean")
hist(distance.matrix)