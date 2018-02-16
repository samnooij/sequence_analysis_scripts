#! /usr/bin/env python

#New Sapovirus metadata parser
from Bio import SeqIO
from sys import argv
import pandas as pd

script, input_file = argv
input_table = pd.read_csv(input_file, sep = '\t')

def parse_genogroup():
    pass

def parse_host():
    pass

def parse_country():
    pass

def parse_year():
    pass

def create_new_table():
    """
    Takes a table with accession IDs in the following format:
    accession_id    genogroup   size_(nt)   strain  host
    
    To create a new table in the format:
    accession_id    genogroup   size_nt strain  host    country year
    
    by reading the table and parsing GenBank.
    
    INPUT = table as pandas dataframe
    OUTPUT = table as pandas dataframe (with more columns)
    """
    for accession in input_table["accession_id"]:
        pass
    


if __name__ == "__main__":
    

OLD_TABLE:



NEW_TABLE:



parse these from:
accession_id = tabular text file
genogroup = tabular text file  (experiment with genbank feature table, "/note=")
size_nt = tabular text file
host = tabular text file (experiment with genbank feature table, "/country=")
country = genbank feature table, "/country="
year = genbank feature table, "/collection_date="

#Note that this would work for HM002617
# It might not work for others, so prepare to throw in some 
try:
# and
if:
# statements
