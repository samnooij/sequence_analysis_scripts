#! /usr/bin/python

#Rename the fasta file to remove complicated names

from Bio import SeqIO

for seq_record in SeqIO.parse("Sapovirus-VP1-translatorx.aaseqs.fasta", "fasta"):
    print ">%s\n%s" % (seq_record.id.split("_")[0], seq_record.seq)
