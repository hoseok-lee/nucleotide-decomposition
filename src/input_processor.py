from sys import stderr
from os.path import splitext, exists

from src.nucleotide_sequence import NucleotideSequence



"""
Given a filename/filepath for a FASTA file, produce a list of custom-defined
nucleotide sequence objects.

Parameters
----------
filename : str
    The filename of the FASTA file

Returns
-------
list[NucleotideSequence]
    List of NucleotideSequence objects

Raises
------
FileError
    Raised when the file is not found
"""
def process_fasta_file (filename: str) -> list[NucleotideSequence]:
    # Check if file exists
    if not exists(filename):
        raise FileNotFoundError

    # Read the entire file in one string
    with open(filename, "r") as file:
        raw = file.read()

    # Filter out comments
    raw = remove_comments(raw)



    # Read the entire file and then split by ">"
    # ">" is the most commonly used delimiter for FASTA files
    # Remove empty strings (and new line characters) in list
    raw = list(filter(
        lambda line: (line and line != '\n'), 
        raw.split('>')
    ))

    # Split each string "chunk" between name and sequence
    # First line will always be the name/description and all following will be
    # the nucleotide (or amino acid) sequence
    raw = list(map(
        lambda fasta_chunk: 
            [
                # Name (first line only)
                fasta_chunk.split('\n')[0], 
                # Sequence (every subsequent line)
                ''.join(fasta_chunk.split('\n')[1:])
            ],
        raw
    ))

    # Generate nucleotide sequence classes
    sequences = list(map(
        lambda fasta_chunk:
            # Name, sequence, and whether the sequence is an amino acid
            # (.faa is the only file extension for amino acid sequences)
            # Amino acid sequences must be converted to nucleotide sequences
            # before the reverse complement is computed
            NucleotideSequence(
                name=fasta_chunk[0], 
                sequence=fasta_chunk[1], 
                amino_acid=(process_file_ext(filename) == '.faa')
            ),
        raw
    ))
    
    return sequences

def remove_comments (file: str) -> str:
    # Filter out all comments i.e. new lines starting with a semicolon ";"
    return '\n'.join(list(filter(
        lambda line: (line.startswith(';') == False), 
        file.split('\n')
    )))

def process_file_ext (filename: str) -> str:
    # Retrieve file extension
    return splitext(filename)[1]