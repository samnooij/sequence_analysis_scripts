#! /usr/bin/env python3

# Author: Sam Nooij
# Institute: Leiden University Medical Center (LUMC)
# Date: 2024-01-25

import argparse
import copy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    """
    parse arguments from the command-line:
        -i/--input = input file (GenBank format)
        -s/--start = gene to start from
        -e/--end = gene to end with
        -o/--output = output file (default = "extracted_records.gbff")
    """
    parser = argparse.ArgumentParser(
        description = "Search a a GenBank file for a start and end gene of interest"
        ", then return only the records within that range."
        "The output is again a GenBank file.")
    
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
        "-s",
        "--start",
        dest = "start",
        required = True,
        type = str,
        help = "Gene (locus tag) to use as starting gene"
        )
    
    required.add_argument(
        "-e",
        "--end",
        dest = "end",
        required = True,
        type = str,
        help = "Gene (locus tag) at the end of the cluster")
    
    optional = parser.add_argument_group("Optional arguments")
    
    optional.add_argument(
        "-o",
        "--output",
        dest = "output",
        required = False,
        type = str,
        default = "extracted_records.gbff",
        help = "Output filename (defualt: extracted_records.gbff)"
        )

    args = parser.parse_args()
    
    return args

def modify_feature(gb_feature, start_pos):
    """
Take a GenBank feature and start position as input,
return the feature with modified (- start position) location.
    """
    updated_feature = copy.deepcopy(gb_feature)
    
    new_start = gb_feature.location.start - start_pos
    new_end = gb_feature.location.end - start_pos
    new_feature_location = FeatureLocation(new_start, new_end, gb_feature.location.strand)
    updated_feature.location = new_feature_location
    
    return updated_feature


def extract_features(inputfile, start_locus, end_locus):
    """
    Parse a GenBank file to extract only the genes of interest, marked by the
    starting and ending locus tag as entered as command-line arguments.
    Return: GenBank record with only the features of interest
    """
    for record in SeqIO.parse(inputfile, "gb"):
    
        sequence_id = record.id
        sequence_name = record.name
        sequence_description = record.description
        original_sequence = record.seq
        
        features_of_interest = []
        
        print("Original sequence length: %s bp" % len(original_sequence))
        
        start_found = False
        
        for feature in record.features:
            if start_found:
                features_of_interest.append(modify_feature(feature, start_position))
            else:
                pass
    
            for key, value in feature.qualifiers.items():
                if key == "locus_tag" and value[0] == start_locus:
                    start_position = int(feature.location.start)
    
                    print("Start is: %s" % start_position)
                    
                    features_of_interest.append(modify_feature(feature, start_position))
                    
                    start_found = True
                
                if start_found and key == "locus_tag" and value[0] == end_locus:
                    end_position = int(feature.location.end)
                    
                    features_of_interest.append(modify_feature(feature, start_position))
                    
                    start_found = False
                    break
    
        try:
            extracted_sequence = original_sequence[start_position:end_position]
            print("Extracted sequence length: %s bp" % len(extracted_sequence))
            
            extracted_record = SeqRecord(extracted_sequence)
            extracted_record.id = sequence_id
            extracted_record.name = sequence_name
            extracted_record.description = sequence_description
            extracted_record.features = features_of_interest
            
            print("Your record of interest is now:")
            print(extracted_record)
            
            return extracted_record

        except UnboundLocalError:
            print("Start gene not found! (Did you make a typo?)")
            break
        
    return None


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
        "    START: {1}\n"
        "    END: {2}\n"
        "    OUTPUT: {3}\n".format(
            arguments.input,
            arguments.start,
            arguments.end,
            arguments.output
            )
        )
    
    print(message)
    
    extracted_record = extract_features(
        inputfile = arguments.input,
        start_locus = arguments.start,
        end_locus = arguments.end
    )
    
    if type(extracted_record) != SeqRecord:
        print("Could not extract the features of interest...")
        return 2

    SeqIO.write(extracted_record, arguments.output, "gb")
    print("Wrote the extracted records to %s !" % arguments.output)
    return 0
    
# EXECUTE script--------------------------------------------
if __name__ == "__main__":
    exit(main())
