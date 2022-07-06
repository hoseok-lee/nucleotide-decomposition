import argparse
from sys import stderr
from os.path import join, dirname

import altair as alt
import pandas as pd
from scipy.spatial.distance import cosine, euclidean

from src.fasta_processor import open_fasta_file, create_fasta_file



if __name__ == "__main__":
    # Read provided arguments
    parser = argparse.ArgumentParser(description='Given a FASTA file of \
        nucleotide (or amino acid) sequences, compute their reverse \
        complements and graph their composition of nucleotide bases (or amino \
        acids). Compare the composition of the first two sequences. \
        Note that in verbose mode, the similarity comparison may be inaccurate \
        due to the nature of Euclidean distances in higher dimensions.')
    parser.add_argument('input', help='an input FASTA file')
    parser.add_argument('output', help='name for the output FASTA file')
    parser.add_argument('-v', 
        help='verbose decomposition (avoid estimation of probabilistic bases)',
        action='store_true'
    )
    args = vars(parser.parse_args())
    input_filename, output_filename, verbose = args.values()



    '''
        REVERSE COMPLEMENT AND DECOMPOSITION
        Open the FASTA file, compute the reverse complement
        Store the reverse complements in the output file
        Compute the compositions of the reverse complements
    '''
    # Open the file and generate the sequences
    input_filename = join(dirname(__file__), input_filename)
    sequences = open_fasta_file(input_filename)

    # Generate the reverse complements
    reverse_complements = list(map(
        lambda sequence: sequence.reverse_complement(),
        sequences
    ))

    # Write reverse complements to new file
    output_filename = join(dirname(__file__), output_filename)
    create_fasta_file(output_filename, reverse_complements)

    # Compute the compositions
    compositions = dict(map(
        lambda sequence: 
            (
                sequence.get_name(),
                sequence.composition(verbose=verbose)
            ),
        reverse_complements
    ))



    '''
        GRAPHICAL REPRESENTATION
        Trellis bar plots with Altair and pandas
    '''
    # Generate pandas Dataframe
    # Transpose to prepare for melt and add column for the name of sequence
    source = pd.DataFrame(compositions).T
    source = source.rename_axis('name').reset_index()
    # Melt the DataFrame and rename the default columns to base and composition
    source = pd.melt(source, id_vars='name').rename(
        columns={'variable': 'base', 'value': 'composition'})
    
    # Use Altair to generate a standard bar graph
    bars = alt.Chart(source).mark_bar().encode(
        x=alt.X('composition:Q', title=None, axis=alt.Axis(format='%')),
        y=alt.Y('base:N', title=None),
        color='base:N'
    )

    # Draw text on each bar as a rounded percentage
    text = bars.mark_text(
        align='left',
        baseline='middle',
        dx=5
    ).encode(
        text=alt.Text('composition:Q', format='.0%')
    )

    # Turn into a Trellis plot with each sequence as a facet
    chart = (bars + text).facet(
        facet=alt.Facet(
            'name:N', 
            title='Nucleotide Sequence Decomposition',
            header=alt.Header(labelLimit=400)
        ),
        columns=2
    )

    # Save as HTML
    chart.save('chart.html')



    '''
        SEQUENCE COMPOSITION SIMILARITY
        Cosine similarity and Euclidean distance as a measurement of difference
    '''
    # Grab the first two compositions and compute the cosine similarity
    assert(len(compositions) >= 2)
    label1, label2 = list(compositions.keys())[:2]
    data1 = list(compositions[label1].values())
    data2 = list(compositions[label2].values())

    assert(len(data1) == len(data2))
    similarity = euclidean(data1, data2)

    # Print result
    print("Euclidean distance of {:.3f} between \
        \n\n\t1. \'{}\' \n\t2. \'{}\'\n".format(similarity, label1, label2), 
        file=stderr)
    print("The two sequences are {}significantly different."\
        .format("" if similarity > 0.1 else "not "), 
        file=stderr)