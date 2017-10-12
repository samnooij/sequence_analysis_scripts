#! /usr/bin/python

#Add genogroup labels to the beginning of the amino acid fasta.

from Bio import SeqIO

directory = "../data/"

nt_file = directory + "Sapovirus_VP1-manually_sliced.fasta"
aa_file = directory + "Sapovirus_VP1-manually_sliced_aa.fasta"

translator_dict = {}

for seq_record in SeqIO.parse(nt_file, "fasta"):
    genogroup = seq_record.id.split('|')[0]
    id = seq_record.id.split('|')[1]
    translator_dict[id] = genogroup

#Make a dictionary with accession IDs as keys and genogroups
# as values to be able to add genogroups to the names in the
# fasta file. E.g. 'AB221130.1': 'GVII'

with open(directory + "Sapovirus_VP1-manually_sliced_aa+genogroups.fasta", 'w') as outfile:
    for seq_record in SeqIO.parse(aa_file, "fasta"):
        old_name = seq_record.description
        accession_id = seq_record.id.split('_')[0]
        #The translated file tells the reading frame (_1) in
        # the name line, which I split out here.
        new_name = "%s|%s" % (translator_dict[accession_id], old_name)
        outfile.write(">%s\n%s\n" % (new_name, seq_record.seq))
