#! /usr/bin/python

#Separate VP1 sequences from different genogroups into 
# different fasta files.

import glob
import os
from Bio import SeqIO

seq_directory = "../data/"

aa_file = seq_directory + "Sapovirus_VP1-manually_sliced_aa+genogroups.fasta"
nt_file = seq_directory + "Sapovirus_VP1-manually_sliced.fasta"

genogroup_dict_aa = {}

for seq_record in SeqIO.parse(aa_file, "fasta"):
    genogroup = seq_record.id.split('|')[0]
    if genogroup not in genogroup_dict_aa:
        genogroup_dict_aa[genogroup] = []
        genogroup_dict_aa[genogroup].append(seq_record)
    else:
        genogroup_dict_aa[genogroup].append(seq_record)
#Collect all sequence records per genogroup in a dictionary.

print("%i Genogroups have been observed in the aa file." % len(genogroup_dict_aa))
print("And these are now written to the following files:")

if not os.path.exists(seq_directory + "separate-aa/"):
    os.mkdir(seq_directory + "separate-aa/")
else:
    pass

for group in genogroup_dict_aa.keys():
    filename = "%sseparate-aa/Sapovirus_VP1_aa-%s.fasta" % (seq_directory, group)
    SeqIO.write(genogroup_dict_aa[group], filename, "fasta")
    print(filename)
#And write these to separate files.        

genogroup_dict_nt = {}

for seq_record in SeqIO.parse(nt_file, "fasta"):
    genogroup = seq_record.id.split('|')[0]
    if genogroup not in genogroup_dict_nt:
       genogroup_dict_nt[genogroup] = []
       genogroup_dict_nt[genogroup].append(seq_record)
    else:
       genogroup_dict_nt[genogroup].append(seq_record)
       

print("%i Genogroups have been observed in the nt file." % len(genogroup_dict_nt))
print("And these are now written to the following files:")

if not os.path.exists(seq_directory + "separate-nt/"):
    os.mkdir(seq_directory + "separate-nt/")
else:
    pass

for group in genogroup_dict_nt:
    filename = "%sseparate-nt/Sapovirus_VP1_nt-%s.fasta" % (seq_directory, group)
    SeqIO.write(genogroup_dict_nt[group], filename, "fasta")
    print(filename)