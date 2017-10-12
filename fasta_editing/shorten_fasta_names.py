#! /usr/bin/python

#Rename the fasta file to remove complicated names

from Bio import SeqIO

file = "Sapovirus-VP1+genogroups.fasta"

for seq_record in SeqIO.parse(file, "fasta"):
    print ">%s\n%s" % (seq_record.id, seq_record.seq)
