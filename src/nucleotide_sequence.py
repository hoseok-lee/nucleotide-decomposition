from __future__ import annotations
from collections import Counter, defaultdict
from functools import reduce

from src.constants import INV_AMINO_CODON, COMPL_BASE, SEQ_REPR, STD_BASES



class NucleotideSequence (object):
    name = ""
    sequence = ""
    amino_acid = False

    def __init__ (self, name: str, sequence: str, amino_acid=False) -> None:
        self.name = name
        self.sequence = sequence
        # Convert all lower-case to upper-case

        self.amino_acid = amino_acid

        # Convert amino acids into nucleotide bases
        if amino_acid:
            self.sequence = ''.join(list(map(
                lambda codon: INV_AMINO_CODON[codon],
                self.sequence
            )))

        # Check if its DNA or RNA
        self.uracil_present = False
        if "U" in self.sequence:
            self.uracil_present = True

    def reverse_complement (self) -> NucleotideSequence:
        return NucleotideSequence(
            name=self.name,
            # The map() computes the complementary sequence, which is then 
            # converted to a list and reversed
            sequence=''.join(list(map(
                lambda base: COMPL_BASE[base],
                self.sequence
            ))[::-1])
        )

    def composition (self, mode='estimate') -> int:
        # Generate a Counter that counts all the occurrences of each base
        # len(self.sequence) will return the same as sum(Counter.values())
        occurrence = Counter(self.sequence)

        # Remove all skips
        del occurrence['-']

        # Replace T with U before composition calculation
        # This will be accounted for at the end with uracil presence
        occurrence['U'] = occurrence['T']
        del occurrence['T']

        # Partition occurrences into standard bases (A, G, C, T, or U) and non-
        # standard bases
        special_bases, standard_bases = reduce(
            lambda part, key_value:
                part[key_value[0] in STD_BASES].update([key_value]) or part,
            occurrence.items(),
            ({}, {})
        )

        # Simply return the occurrence if not required to estimate
        if mode == 'standard':
            return occurrence

        # Cut all special bases
        if mode == 'truncate':
            return standard_bases



        # ESTIMATED COMPOSITION
        # Special bases, such as N (which can be any base) or S (which can only
        # be C or G) for example, can be used to estimate the composition by 
        # distributing the occurrence along the possible constituents
        # For example, if the sequence contains 5 S bases, then there's a 
        # possibility of a mixture of 5 C or G bases, leading to an estimated
        # composition of 2.5 C and 2.5 G
        estimated_comp = dict(reduce(
            lambda comp, key_value:
                reduce(
                    lambda pseudo, base:
                        pseudo.update({
                            base: comp[base] + (
                                key_value[1] / len(SEQ_REPR[key_value[0]])
                                if base in SEQ_REPR[key_value[0]] else 0
                            )
                        }) or pseudo,
                    SEQ_REPR[key_value[0]],
                    defaultdict(int)
                ),
            special_bases.items(),
            standard_bases
        ))

        # Replace U with T if DNA
        if not self.uracil_present:
            estimated_comp['T'] = estimated_comp['U']
            del estimated_comp['U']

        return estimated_comp

    def get_name (self) -> str:
        return self.name

    def __repr__ (self) -> str:
        return "<name=\"{}\", sequence={}>".format(self.name, self.sequence)