from nucleotide_sequence import NucleotideSequence

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
KeyError
    when a key error
OtherError
    when an other error
"""
def process_fasta_file (filename: str) -> list[NucleotideSequence]:
    # Read the entire file in one string
    file = open(filename, "r")
    raw = file.read()
    file.close()
    
    # Filter out all comments i.e. lines starting with a semicolon ";"
    raw = '\n'.join(list(filter(
        lambda line: (line.startswith(';') == False), 
        raw.split('\n')
    )))

    # Read the entire file and then split by ">"
    # ">" is the most commonly used delimiter for FASTA files
    # Remove empty strings (and new line characters) in list
    raw = list(filter(
        lambda sequence: (sequence and sequence != '\n'), 
        raw.split('>')
    ))

    # Generate nucleotide sequence classes
    sequences = list(map(lambda sequence: NucleotideSequence(sequence), raw))
    return sequences