#! /usr/bin/python

#A script that downloads protein sequences from
# the NCBI that go with the accession IDs in the given file.
# It extracts CDSs from GenBank records to parse the protein
# sequences.
#Note: this script is made for Sapovirus sequences.

#Example use: $ python download_genbank_cds_proteins.py VP1.lst

from sys import argv
from Bio import Entrez
import pandas as pd
import textwrap

Entrez.email = "sam.nooij@rivm.nl" #Tell NCBI who you are!

# Define required files: input and output
filename = argv[1] #first command line argument
output_file = "Sapovirus-%s-proteins.fasta" % filename.split('.')[0]

# Read the accession numbers from the input file
accession_df = pd.read_csv(filename, header=None)
accession_ids = ' '.join([id for id in accession_df[0]])

# Download the associated GenBank records
handle = Entrez.efetch(db="nuccore", id=accession_ids, 
                       rettype = "gb", retmode = "XML")
records = Entrez.read(handle)

# Parse the CDS features from the GenBank records and write
# all protein sequences to a new file
with open(output_file, 'w') as fasta:
    for record in records:
        name = record["GBSeq_primary-accession"]
        protein_sequence = ''
        for cds in records[0]["GBSeq_feature-table"][1:]:
            translation = cds["GBFeature_quals"][-1]["GBQualifier_value"]
            protein_name = cds["GBFeature_quals"][2]["GBQualifier_value"]
            protein_id = cds["GBFeature_quals"][3]["GBQualifier_value"]
            name += "-%s[%s]" % (protein_name, protein_id)
            protein_sequence += translation
        fasta.write(">%s\n" % name)
        fasta.write(textwrap.fill("%s" % protein_sequence, width = 60))
        fasta.write('\n')
