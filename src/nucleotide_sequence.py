from __future__ import annotations

from src.constants import INVERSE_AMINO_CODON, COMPLEMENTARY_BASE




class NucleotideSequence (object):
    name = ""
    sequence = ""
    amino_acid = False

    def __init__ (self, name: str, sequence: str, amino_acid=False) -> None:
        self.name = name
        self.sequence = sequence
        self.amino_acid = amino_acid

        # Convert amino acids into nucleotide bases
        if amino_acid:
            self.sequence = ''.join(list(map(
                lambda codon: INVERSE_AMINO_CODON[codon],
                self.sequence
            )))

    def reverse_complement (self) -> NucleotideSequence:
        return NucleotideSequence(
            name=self.name,
            # The map() computs the complementary sequence, which is then 
            # converted to a list and reversed
            sequence=''.join(list(map(
                lambda base: COMPLEMENTARY_BASE[base],
                self.sequence
            ))[::-1])
        )

    def __repr__ (self) -> str:
        return "<name=\"{}\", sequence={}>".format(self.name, self.sequence)