#! /usr/bin/python

''' 
A simple script that allows the user to produce a shorter version
of GenBank Feature tables. It takes two command line arguments:
the file you want to read and the number of record to write.

EXAMPLE USE:
python first_X_lines.py viral.1.genomic.gbff 100

Author: S. Nooij
Email: sam.nooij [at] rivm [dot] nl
creation date: 27-11-2015
'''

from sys import argv

SCRIPT, INPUTFILE, RECORDS = argv

RECORDS = int(RECORDS)
# RECORDS needs to be a number (integer), not a string

OUTPUTFILE = '.'.join(INPUTFILE.split('.')[:-1])
OUTPUTFILE += "-first%dRecords." % RECORDS
OUTPUTFILE += INPUTFILE.split('.')[-1]

print OUTPUTFILE 

#These lines take care of the naming of the outputfile
# by inserting "-firstXRecords" before the file extension,
# where X is the number of records you specified.

with open(INPUTFILE, 'r') as readfile:
    with open(OUTPUTFILE, 'w') as writefile:
        for line in readfile:
            writefile.write(line)
            if line.startswith('//'):
                RECORDS -= 1
            if RECORDS < 1:
                break

