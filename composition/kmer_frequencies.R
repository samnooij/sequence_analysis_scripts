# Calculate monomer, dimer and trimer frequencies of
# the sequences in a given fasta file.

# Author: Sam Nooij
# Date: 12 Oct 2017
# Email: sam.nooij@rivm.nl

library('seqinr')

sequences <- read.fasta(file = "sequences.fasta")

monomers <- function(sequence) {
  frequencies <- count(seq = sequence, wordsize = 1)
  return(frequencies)
}

dimers <- function(sequence) {
  frequencies <- count(seq = sequence, wordsize = 2)
  return(frequencies)
}

trimers <- function(sequence) {
  frequencies <- count(seq = sequence, wordsize = 3)
  return(frequencies)
}

monomer.df <- as.data.frame(do.call(rbind, lapply(X = sequences, FUN = monomers)))
monomer.df$name <- rownames(monomer.df)
#Reorder to put the 'name' column first:
monomer.df <- monomer.df[c(ncol(monomer.df), 1:ncol(monomer.df) - 1)]

dimer.df <- as.data.frame(do.call(rbind, lapply(X = sequences, FUN = dimers)))
dimer.df$name <- rownames(dimer.df)
#Reorder to put the 'name' column first:
dimer.df <- dimer.df[c(ncol(dimer.df), 1:ncol(dimer.df) - 1)]

trimer.df <- as.data.frame(do.call(rbind, lapply(X = sequences, FUN = trimers)))
trimer.df$name <- rownames(trimer.df)
#Reorder to put the 'name' column first:
trimer.df <- trimer.df[c(ncol(trimer.df), 1:ncol(trimer.df) - 1)]

#Write to .csv files
write.table(x = monomer.df, file = "monomer_frequencies.csv", sep = ",", row.names = FALSE, col.names = TRUE)
write.table(x = dimer.df, file = "dimer_frequencies.csv", sep = ",", row.names = FALSE, col.names = TRUE)
write.table(x = trimer.df, file = "trimer_frequencies.csv", sep = ",", row.names = FALSE, col.names = TRUE)
