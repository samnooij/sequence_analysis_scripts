#! /usr/bin/python

#Date: 04-03-2016
#Author: S.Nooij
#Email: sam.nooij [at] rivm [dot] nl

from oldowan.genbank import iterate_genbank
from Bio import SeqIO
from Bio.Seq import Seq #for testing translations
from Bio.SeqRecord import SeqRecord

INPUTFILE = "Reference_sequences.gb"
KEYWORDS = ["hn", "hemagglutinin", "neuraminidase"]
SEQ_OUTPUT = "tmp-Mumps_HN_genes_from_references.fasta"
TRANS_OUTPUT = "tmp-Mumps_HN_translations_from_references.fasta"

def slice_gene(seq, pos):
    """
    Slices the sequence of the gene from
    the complete sequence.
    seqe = sequence (string or Seq object)
    pos = list of begin and end positions
    """
    return seq[pos[0]:pos[1] + 1]
    #the + 1 was added, because slicing in Python does not
    # include the "end position"

def make_integers(str_list):
    """
    Make a list of integers from
    a list of strings (position)
    and convert to Python indices
    (- 1).
    """
    new_list = []
    for number in str_list:
        if '>' in number:
            number = number.replace('>', '')
        new_list.append(int(number)-1)
    return new_list

    #Do not use for HN - "complete CDS" sequences do not
    # include UTRs.
def add_utrs(pos_list):
    """
    Add the UTRs to the positions of the SH sequence.
    INPUT: list of 2 ints
    OUTPUT: list of 2 ints (changed)
    """
    five_prime = 50
    three_prime = 92
    new_positions = [pos_list[0]-five_prime, \
    pos_list[1]+three_prime]
    return new_positions

if __name__ == '__main__':
    hn_sequences = []
    all_accessions = []
    hn_accessions = []
    hn_translations = []

    for entry in iterate_genbank(INPUTFILE):
        hn = False
        accession = entry['accession']
        all_accessions.append(accession)
        features = entry['features']
        if features:
            source = features[0]
            for element in features:
                if element[0] == 'CDS':
                    cds = element
                    cds_content = cds[1]
                    position = make_integers(cds_content.pop(0).split('..'))
                    product = dict(cds_content).get('product')
                    translation = Seq(dict(cds_content).get('translation').strip('"'))
                    for words in KEYWORDS:
                        if words in product.lower():
                            hn = True
#                        print "Pos: %s,\tProduct: %s" % (position, product)
                if hn:
                    break

        sequence = Seq(entry['sequence'])

        if hn:
            hn_accessions.append(accession)
            hn_sequence = slice_gene(sequence, position)
            product = product.replace(' ','_')
### for debugging
#            print accession, position, product, translation[:10], sequence[:10]
#            print hn_sequence[:10], translation[-10:], hn_sequence.translate()[-10:]

            seq_record = SeqRecord(hn_sequence)
            seq_record.id = "%s_%s" % (accession, product)
            seq_record.description = "Mumps virus, nt sequence"
            trans_record = SeqRecord(translation)
            trans_record.id = "%s_%s" % (accession, product)
            trans_record.description = "Mumps virus, AA sequence"
            hn_sequences.append(seq_record)
            hn_translations.append(trans_record)
            print accession, len(hn_sequence)
        hn = False
        
        if accession == "AB000388":
            hn_accessions.append(accession)
            position = [6614, 8362]
            #Note that these positions correspond to those
            # listed in the Excel sheet
            hn_sequence = slice_gene(sequence, position)
            product = "HN"
            seq_record = SeqRecord(hn_sequence)
            seq_record.id = "%s_%s" % (accession, product)
            seq_record.description = "Mumps virus, nt sequence"
            hn_sequences.append(seq_record)
            
    with open(SEQ_OUTPUT, 'w') as seq_file:
        SeqIO.write(hn_sequences, seq_file, "fasta")
    with open(TRANS_OUTPUT, 'w') as trans_file:
        SeqIO.write(hn_translations, trans_file, "fasta")
