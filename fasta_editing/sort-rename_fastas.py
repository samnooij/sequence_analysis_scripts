#! /usr/bin/python

#6-10-2017
#Sort Sapovirus sequences according to genogroup number

from Bio import SeqIO
import textwrap

nt_file = "data/Sapovirus_VP1-manually_sliced.fasta"
aa_file = "data/Sapovirus_VP1_aa.fasta"

nt_file_new = "data/Sapovirus_VP1_nt_sorted.fasta"
aa_file_new = "data/Sapovirus_VP1_aa_sorted.fasta"


nt_list = []
aa_list = []

for seq_record in SeqIO.parse(nt_file, "fasta"):
    nt_list.append((seq_record.id, seq_record.seq))

for seq_record in SeqIO.parse(aa_file, "fasta"):
    aa_list.append((seq_record.id, seq_record.seq))

sort_order = {'GI' : 1, 'GII' : 2, 'GIII' : 3, 'GIV' : 4, 'GV' : 5, 'GVI' : 6,
              'GVII' : 7, 'GVIII' : 8, 'GIX' : 9, 'GX' : 10, 'GXI' : 11, 
              'GXII' : 12, 'GXIII' : 13, 'GXIV' : 14, 'GXV' : 15}

nt_list.sort(key=lambda val: sort_order[val[0].split("|")[0]])
aa_list.sort(key=lambda val: sort_order[val[0].split("|")[0]])

with open(nt_file_new, 'w') as new_nt:
    for entry in nt_list:
        new_nt.write(">%s\n" % entry[0])
        new_nt.write(textwrap.fill("%s" % entry[1], width = 60))
        new_nt.write("\n")
        
with open(aa_file_new, 'w') as new_aa:
    for entry in aa_list:
        new_aa.write(">%s\n" % entry[0])
        new_aa.write(textwrap.fill("%s\n" % entry[1], width = 60))
        new_aa.write("\n")