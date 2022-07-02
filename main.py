from audioop import reverse
import sys
from os.path import join, dirname

from src.input_processor import process_fasta_file
from src.nucleotide_sequence import NucleotideSequence



if __name__ == "__main__":
    # Main control loop
    filename = join(dirname(__file__), "nucleotide-sample.fasta")
    sequences = process_fasta_file(filename)
    reverse_complements = list(map(
        lambda sequence: sequence.reverse_complement(),
        sequences
    ))

    