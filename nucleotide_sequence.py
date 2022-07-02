class NucleotideSequence (object):
    name = ""
    sequence = ""

    def __init__(self, raw_fasta: str):
        # Split the raw FASTA string between name and sequence
        # First line will always be the name/description and all following
        # will be the nucleotide sequence
        raw_fasta = raw_fasta.split('\n')
        self.name = raw_fasta[0]
        self.sequence = ''.join(raw_fasta[1:])