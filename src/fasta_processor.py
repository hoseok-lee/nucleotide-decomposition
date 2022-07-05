from os import remove
from sys import stderr
from os.path import splitext, exists

from src.nucleotide_sequence import NucleotideSequence



def open_fasta_file (filename: str) -> list[NucleotideSequence]:
    # Read the entire file in one string
    with open(filename, "r") as file:
        raw = file.read()

    # Decode FASTA chunks
    fasta_chunks = decode_fasta_chunks(raw)

    # Generate nucleotide sequence classes
    sequences = list(map(
        lambda key_value:
            # Name, sequence, and whether the sequence is an amino acid
            # (.faa is the only file extension for amino acid sequences)
            # Amino acid sequences must be converted to nucleotide sequences
            # before the reverse complement is computed
            NucleotideSequence(
                name=key_value[0], 
                sequence=key_value[1],
                amino_acid=(get_file_ext(filename) == '.faa'),
                uracil_present=('U' in key_value[1])
            ),
        fasta_chunks.items()
    ))
    
    return sequences

def create_fasta_file (filename: str, sequences: list[NucleotideSequence]) -> None:
    # Create dictionary of names and sequence strings
    fasta_chunks = dict(map(
        lambda sequence: (sequence.get_name(), sequence.get_sequence()),
        sequences
    ))

    # Encode FASTA chunks
    raw = encode_fasta_chunks(fasta_chunks)

    # Open file for writing
    with open(filename, "w") as file:
        file.write(raw)



def decode_fasta_chunks (file: str) -> dict[str, str]:
    # Filter out all comments i.e. new lines starting with a semicolon ";"
    file = '\n'.join(list(filter(
        lambda line: (line.startswith(';') == False), 
        file.split('\n')
    )))

    # Read the entire file and then split by ">"
    # ">" is the most commonly used delimiter for FASTA files
    # Remove empty strings (and new line characters) in list
    file = list(filter(
        lambda line: (line and line != '\n'), 
        file.split('>')
    ))

    # Split each string "chunk" between name and sequence
    # First line will always be the name/description and all following will be
    # the nucleotide (or amino acid) sequence
    fasta_chunks = dict(map(
        lambda fasta_chunk: 
            (
                # Name (first line only)
                fasta_chunk.split('\n')[0], 
                # Sequence (every subsequent line)
                ''.join(fasta_chunk.split('\n')[1:])
            ),
        file
    ))

    return fasta_chunks

def encode_fasta_chunks (fasta_chunks: dict[str, str]) -> str:
    # Join all FASTA chunks by splitting the name and sequence with a new-line
    # and initializing with the delimiter ">"
    # Note that the sequence can be split over any number of lines (easiest
    # would be on one single line), and its end is only defined by the start of
    # the next sequence with ">"
    file = '\n'.join(list(map(
        lambda key_value: ">{}\n{}\n".format(key_value[0], key_value[1]),
        fasta_chunks.items()
    )))

    return file

def get_file_ext (filename: str) -> str:
    # Retrieve file extension
    file_ext = splitext(filename)[1]
    return file_ext