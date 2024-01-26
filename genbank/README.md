# GenBank parser scripts (mostly using Biopython)

Scripts that I have used to either edit/parse GenBank files or download files from GenBank.

## GenBank file parsing scripts

(Both of these use [Biopython](https://biopython.org/) - so make sure that library is available before using!)

  - [genbank_extractor](genbank/genbank_extractor.py): This script reads a file in GenBank format (for example, a genome annotation file - like those produced by [Bakta](https://bakta.readthedocs.io/en/latest/) for prokaryotic genome sequences), and takes as extra input arguments a start and end gene to extract only the features from the start to the end. This has been used to extract gene clusters of interest for visualisation with [clinker](https://github.com/gamcil/clinker). The script uses the argparse library to handle command-line arguments and provide help with the `-h` flag.

  - [genbank_to_reoriented_fasta](genbank/genbank_to_reoriented_fasta.py): This script takes a GenBank file (again, a genome annotation file) as input together with a gene locus tag to then extract the sequence and re-orient it to start at the provided gene. This can then be used as input to Bakta to produce a annotation file that starts at the desired gene. (It should be possible to sort the features in the GenBank file and change their positions, but I could not figure out how to do that so easily, so I settled for this workaroud.)
