import sys
from os.path import join, dirname
import altair as alt
from vega_datasets import data

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


    compositions = dict(map(
        lambda sequence: 
            (
                sequence.get_name(),
                sequence.composition()
            ),
        reverse_complements
    ))

    print(compositions)

    #source = data.cars()
    #alt.Chart(source).mark_bar().encode(
    #    alt.X("Horsepower:Q", bin=True),
    #    y='count()',
    #    row='Origin'
    #)