#! /usr/bin/python

#Author: Sam Nooij
#Email: sam.nooij@rivm.nl
#Date: 28-11-2016

from __future__ import print_function
import urllib
import numpy as np
import glob
import os
import pandas as pd
from Bio import SeqIO
from Bio import Entrez


VIRUS_OF_INTEREST = "mumps virus"
EMAIL_ADDRESS = "sam.nooij@rivm.nl"
#These are now hard-coded, but maybe it is nicer to
# move them to a config file or something.

Entrez.email = EMAIL_ADDRESS

print("Starting the collection of sequences of your virus of interest from the NCBI.")
print("Your virus is: %s" % VIRUS_OF_INTEREST)
print("First the TaxID is looked up.\n")

## STEP 1: retrieve the TaxID
def retrieve_taxid(virus_name):
    """Given the name of a virus, returns the TaxID from Entrez."""
    handle = Entrez.esearch(db="taxonomy", term=VIRUS_OF_INTEREST)
    record = Entrez.read(handle)
    taxid = record["IdList"][0]
    handle.close()
	
    return taxid

taxid = retrieve_taxid(VIRUS_OF_INTEREST)

print("The TaxID is: %s" % taxid)

print("\nThen a list of IDs is collected for this TaxID.")

## STEP 2: download sequences belonging to that TaxID
def collect_nucleotide_ids_for_taxid(taxid):
    """Given a TaxID, returns all the nucleotide IDs from the NCBI."""
    handle = Entrez.esearch(db="nucleotide", term="txid%s[Organism]" % taxid, retmax="1000000")
    record = Entrez.read(handle)
    recordlist = record["IdList"]
    length = len(recordlist)
    handle.close()
    print("Number of sequences available: " + str(length))

    return recordlist

recordlist = collect_nucleotide_ids_for_taxid(taxid)

print("These are downloaded in batches of about 200 sequences each.")

def download_nt_sequences_for_ids(
        nucleotide_ids, file_prefix):
    """Given a list of NCBI nt IDs, downloads sequences in batches of 200 sequences each."""
#eutils accepts up to 200 IDs at a time,
# so you may chop the list up to download
# more sequences in batches.
    batch_size = 200
    length = len(nucleotide_ids)
    batch_number = length / batch_size
    print("Number of batches to download: %i" % batch_number)

    batched_list = np.array_split(recordlist, batch_number)

    batch_count = 0
    for batch in batched_list:
        batch_count += 1
        batch_string = ','.join(batch)
        urllib.urlretrieve("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\
?db=nucleotide&id=%s&rettype=fasta&retmode=text" % batch_string,
"%s%s.fasta" % (file_prefix, batch_count))
    print("Done downloading batches!")

    return None

sequence_file_prefix = VIRUS_OF_INTEREST.replace(' ', '_') + '_sequences'
download_nt_sequences_for_ids(recordlist, sequence_file_prefix)

def concatenate_fastas(file_prefix):
    """Concatenate fasta files with the same prefix, followed by a number.
E.g. 'mumps_sequences1.fasta' up to 'mumps_sequences100.fasta'.
Also deletes the partial files after concatenating into one file."""
    #concatenate the different fasta files (batches)
    batch_files = glob.glob(file_prefix + "[0-9]*.fasta")

    with open(file_prefix + ".fasta", 'w') as concatenated_fasta:
        for file in batch_files:
            with open(file, 'rU') as part_file:
                seq_records = SeqIO.parse(part_file, 'fasta')
                SeqIO.write(seq_records, concatenated_fasta, 'fasta')

    print("The batches have now been concatenated into one big file.")

    #show the number of sequences in the fasta file
    records = list(SeqIO.parse(file_prefix + ".fasta", "fasta"))
    print("\nNumber of sequences in the concatenated file: %s" % len(records))

    #and delete the separate files
    for	file in batch_files:
        os.remove(file)

    return None

concatenate_fastas(sequence_file_prefix)

def check_fasta(filename):
    """Given a fasta file, parses the first sequence using Biopython's SeqIO."""
    print("\nThese are the details of the first sequence record:")
    for seq_record in SeqIO.parse(filename, "fasta"):
        print("ID = %s" % seq_record.id)
        print("Name = %s" % seq_record.description)
        print("Sequence = %s" % repr(seq_record.seq))
	    #print("Whole record:", seq_record)
        break
        
    return None
    
check_fasta(sequence_file_prefix + '.fasta')

## STEP 3: download the RefSeq reference sequence
#NOTE: this assumes there is only 1 reference sequence,
# which is not always the case...
def download_refseq_reference(taxid):
    """Given a TaxID, downloads the reference sequence from RefSeq to a fasta file."""
    handle = Entrez.esearch(db="nucleotide", 
    term="txid%s[Organism] AND refseq[filter]" %taxid)
    #alternatively, the following query can be used:
    #"%s[ORGN] AND srcdb_refseq[Properties]" % VIRUS_OF_INTEREST
    record = Entrez.read(handle)
    refseq_id = record["IdList"]
    handle.close()
    
    fetched = Entrez.efetch(db="nucleotide", id=refseq_id, 
    rettype="fasta", retmode="text")
    refseq_sequence = fetched.read()
    output_file = "%s-reference.fasta" % VIRUS_OF_INTEREST\
    .replace(' ','_')
    
    with open(output_file,'w') as write_me:
        write_me.write(refseq_sequence)
        
    print("\nThe reference sequence has been downloaded to: \
%s" % output_file)
    
    for seq_record in SeqIO.parse(output_file, 'fasta'):
        print("Name of the reference: %s" 
        % seq_record.description)

    return None

download_refseq_reference(taxid)

## STEP 4: collect the WHO references
#(The accession IDs of the reference sequences
# have been manually written in a textfile.)
#source = http://www.who.int/wer/2012/wer8722.pdf?ua=1

## THESE LINES NEED GENERALIZATION   
directory = "../data/"
reference_file = directory + "Reference_GenBank_IDs.txt"

def read_reference_genotypes(reference_table_file):
    """Reads a table of accession IDs and genotype names.
    Returns a dictionary of these IDs and names."""
    reference_table = pd.read_csv(reference_table_file, delimiter='\t', header=None)

    reference_table.rename(columns={0: 'accession', 1: 'genotype', 2: 'genes', 3: 'strain'}, inplace=True)

    print("Number of sequences: %i" % len(reference_table.index))

    print(reference_table.head())

    reference_id_dict = {}
    for index, row in reference_table.iterrows():
        accession = row[0]
        genotype = row[1]
        #print "Accession: %s\tgenotype: %s" % (accession, genotype)
        reference_id_dict[accession] = 'MuV-' + genotype
 
    #contains_sh = reference_table[reference_table['genes'].str.contains('SH')]
    #contains_hn = reference_table[reference_table['genes'].str.contains('HN')]
    
    return reference_id_dict

reference_ids_dictionary = read_reference_genotypes(reference_file)

def add_genotype_names_to_fasta(
            file_prefix,
            reference_dictionary):
    """Adds the genotype name to the corresponding \
    sequences in the fasta file with all sequences."""
    print("\nThese reference genotypes are now added to \
the fasta with all sequences. The names are added to \
the beginning of the sequence name.")
    sequences_found = 0
    new_fasta = []
    with open(file_prefix + "+genotypes.fasta", 'w') as new_file:
        for seq_record in SeqIO.parse(file_prefix + \
        ".fasta", "fasta"):
            if seq_record.id.split('.')[0] in \
            reference_dictionary:
                sequences_found += 1
                genotype_to_add = \
                reference_dictionary[seq_record.id.split('.')[0]]
                seq_record.id = genotype_to_add + '|' + \
                seq_record.id
                seq_record.description = genotype_to_add + \
                '|' + seq_record.description
            new_file.write('>%s\n%s\n' % \
            (seq_record.description, seq_record.seq))
    
    print("\n==========")
    print("%i reference sequences found and renamed." % 
    sequences_found)
    print("==========")
    
    return None
    
add_genotype_names_to_fasta(sequence_file_prefix, 
reference_ids_dictionary)

#if __name__ == '__main__':   