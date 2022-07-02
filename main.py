import sys

from input_processor import process_fasta_file
from nucleotide_sequence import NucleotideSequence


if __name__ == "__main__":
    # Main control loop
    process_fasta_file("nucleotide-sample.fasta")