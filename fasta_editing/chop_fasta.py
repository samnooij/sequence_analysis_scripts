#! /usr/bin/python

''' 
A simple script that allows the user to produce a shorter version
of fasta files. It takes two command line arguments:
the file you want to read and the number of records to write.

EXAMPLE USE:
python chop_fasta.py virus_sequences.fasta 100

Author: S. Nooij
Email: sam.nooij@rivm.nl
creation date: 27-11-2015
'''

from sys import argv
from Bio import SeqIO

SCRIPT, INPUTFILE, RECORDS = argv

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

#From http://biopython.org/wiki/Split_large_file

RECORDS = int(RECORDS)
# RECORDS needs to be a number (integer), not a string

OUTPUTFILE = '.'.join(INPUTFILE.split('.')[:-1])
extension = INPUTFILE.split('.')[-1]

record_iter = SeqIO.parse(open(INPUTFILE), extension)
for i, batch in enumerate(batch_iterator(record_iter, RECORDS)):
    filename = OUTPUTFILE + "_%i.%s" % ((i + 1), extension)
    handle = open(filename, 'w')
    count = SeqIO.write(batch, handle, extension)
    handle.close()
    print("Wrote %i records to %s" % (count, filename))
