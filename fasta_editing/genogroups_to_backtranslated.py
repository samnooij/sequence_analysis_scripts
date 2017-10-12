#! /usr/bin/python

# EMBOSS tranalign removes genogroup annotation from the 
# fasta files. This script will put them back.

import pandas as pd
from Bio import SeqIO
import textwrap
import glob

reference_file = "data/VP1_GenBank_IDs.tsv"
reference_table = pd.read_csv(reference_file, delimiter = '\t')

reference_id_dict = {}
for index, row in reference_table.iterrows():
    accession = row['accession_id']
    genogroup = row['genogroup']
    reference_id_dict[accession] = genogroup

backtranslated_files = glob.glob("tmp/*_backtranslated*.fa") + glob.glob("tmp/*_backtranslated*.fas")

for file in backtranslated_files:
    # add a '-c' to the name to indicate 'corrected' files
    with open("%s-c" % file, 'w') as new_file:
        for seq_record in SeqIO.parse(file, "fasta"):
            if seq_record.id.split('.')[0] in reference_id_dict:
                genogroup_to_add = reference_id_dict[seq_record.id.split('.')[0]]
            print("%s is genogroup %s" % (seq_record.id, genogroup_to_add))
            seq_record.id = genogroup_to_add + '|' + seq_record.id
            seq_record.description = genogroup_to_add + '|' + seq_record.description
            new_file.write(">%s\n" % seq_record.id)
            new_file.write(textwrap.fill("%s" % seq_record.seq, width = 60))
            new_file.write("\n")
