#! /usr/bin/python

# A script to sort fasta files like the reference file.
## ABORT DEVELOPMENT: it appears that sorting is not necessary now.

### CAUTION: UNFINISHED SCRIPT ###

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob

reference_file = "data/Sapovirus_VP1_nt_sorted.fasta"
order_dict = {}
order_number = 0

for seq_record in SeqIO.parse(reference_file, 'fasta'):
    order_dict[seq_record.id] = order_number
    order_number += 1

print "Reference dictionary:\n", order_dict

files_to_sort = glob.glob("tmp/*-c")

print "\n Files to be sorted:\n", files_to_sort

for unordered_file in files_to_sort:
    print "\nThis file is listed:\n", unordered_file, "\n"
    file_to_sort = []
    for seq_record in SeqIO.parse(unordered_file, "fasta"):
        if not '(' in seq_record.id:
            pass
        else:
            seq_record.id = seq_record.id[:seq_record.id.index('(')]
        file_to_sort.append(seq_record)
        
    #file_to_sort.sort(key=lambda record: order_dict[record.id][1])
    
print "Example sorted file:\n", file_to_sort

print "Example ID:\n", file_to_sort[0].id
