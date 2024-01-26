#! /usr/bin/env python3

# Author: Sam Nooij
# Institute: Leiden University Medical Center (LUMC)
# Date: 2024-01-23

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    """
    parse arguments from the command-line:
        -i/--input = input file (GenBank format)
        -g/--gene = gene to start from
        -o/--output = output file (default = "sorted_sequence.fasta")
    """
    parser = argparse.ArgumentParser(
        description = "Search a a GenBank file for a gene of interest within a longer sequence"
        ", then return the sequence reoriented starting at this gene of interest."
        "If the gene is in reverse complement, return the sequence as reverse complement."
        "The output is saved as fasta file.")
    
    required = parser.add_argument_group("Required arguments")
    
    required.add_argument(
        "-i",
        "--input",
        dest = "input",
        required = True,
        type = str,
        help = "Input file to use (GenBank file)"
        )
    
    required.add_argument(
        "-g",
        "--gene",
        dest = "gene",
        required = True,
        type = str,
        help = "Gene to use as starting gene"
        )
    
    optional = parser.add_argument_group("Optional arguments")
    
    optional.add_argument(
        "-o",
        "--output",
        dest = "output",
        required = False,
        type = str,
        default = "sorted_sequence.fasta",
        help = "Output filename (defualt: sorted_sequence.fasta)"
        )

    args = parser.parse_args()
    
    return args


def find_position(inputfile, gene):
    """
    Parse a GenBank file to find the position and orientation of the gene of interest.
    (written as '[start:end](orientation)', e.g. '[1:500](+)' )
    Return: the start position and orientation as + or -
    """
    for record in SeqIO.parse(inputfile, "gb"):

        for feature in record.features:
            for key, value in feature.qualifiers.items():
                if key == "locus_tag" and value[0] == gene:
                    location = feature.location
                    orientation = str(location)[-2]
                    # could also have been extracted using `feature.location.strand`
                    # as 1 or -1 (?)
                    start_position = int(str(location).split(":")[0].strip("["))
                    # this is a basic Python way, but BioPython also has its own
                    # way: `feature.location.start` !

                    print("Start is: %s, orientation is: %s-strand" % (start_position, orientation))
                break # We need this only once, break the loop
    
    return start_position, orientation


def reorient_sequence(inputfile, start, orientation):
    """
    Read a GenBank file, convert it to a sequence and reorient based on the
    input start position and orientation (reverse complement if orientation
    is --strand (negative strand).
    Return: a SeqRecord object.
    """
    for record in SeqIO.parse(inputfile, "gb"):

        sequence_id = record.id
        sequence_name = record.name
        sequence_description = record.description
        original_sequence = record.seq

        print("Original sequence length: %s bp" % len(original_sequence))

        flipped_sequence = original_sequence[start:] + original_sequence[:start]

        print("Flipped sequence length: %s bp" % len(flipped_sequence))

        if orientation == "-":
            flipped_sequence = flipped_sequence.reverse_complement()
        else:
            pass
            
        flipped_record = SeqRecord(flipped_sequence)
        flipped_record.id = sequence_id
        flipped_record.name = sequence_name
        flipped_record.description = sequence_description
    
    return flipped_record


def main():
    """
    Main execution of the script
    """
    # 1. parse command-line arguments
    arguments = parse_arguments()
    
    message = (
        "\n"
        "These are the input arguments you have provided:\n"
        "    INPUT: {0}\n"
        "    GENE: {1}\n"
        "    OUTPUT: {2}\n".format(
            arguments.input,
            arguments.gene,
            arguments.output
            )
        )
    
    print(message)
    
    start_position, orientation = find_position(inputfile = arguments.input,
                                                gene = arguments.gene)
    
    flipped_record = reorient_sequence(inputfile = arguments.input,
                                       start = start_position,
                                       orientation = orientation)
    
    SeqIO.write(flipped_record, arguments.output, "fasta")
    print("Wrote the reoriented sequence to %s !" % arguments.output)
    
# EXECUTE script--------------------------------------------
if __name__ == "__main__":
    main()
