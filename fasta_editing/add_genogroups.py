#! /usr/bin/python

#A script to add genogroup annotations to a fasta file,
# using a hand-made annotation file and the accession IDs
# in the fasta file.
#Example use: $ python add_genogroups.py VP1

#Taken from the 4-1-2017 folder, modified for use in the current directory

import pandas as pd
from Bio import SeqIO
import textwrap

reference_file = "data/VP1_GenBank_IDs.tsv"
reference_table = pd.read_csv(reference_file, delimiter='\t')

print("Number of sequences: %i" % len(reference_table.index))

reference_id_dict = {}
for index, row in reference_table.iterrows():
    accession = row['accession_id']
    genogroup = row['genogroup']
    reference_id_dict[accession] = genogroup

with open("data/Sapovirus_VP1_aa.fasta",
 'w') as new_file:
    for seq_record in SeqIO.parse("data/Sapovirus_VP1-manually_sliced_aa.fasta", "fasta"):
        if seq_record.id.split('.')[0] in reference_id_dict:
            genogroup_to_add = reference_id_dict[seq_record.id.split('.')[0]]
            print("%s is genogroup %s" % (seq_record.id, genogroup_to_add))
            seq_record.id = genogroup_to_add + '|' + seq_record.id
            seq_record.description = genogroup_to_add + '|' + seq_record.description
        new_file.write(">%s\n" % seq_record.id)
        new_file.write(textwrap.fill("%s" % seq_record.seq, width = 60))
        new_file.write("\n")
