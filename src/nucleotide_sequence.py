from __future__ import annotations
from collections import Counter
from functools import reduce

from src.constants import INV_AMINO_CODON, COMPL_BASE, SEQ_REPR, STD_BASES



"""
A class representing a single nucleotide sequence. It stores the name and 
sequence (string of nucleotide bases), performing operations such as conversion
from amino acid sequences to nucleotide sequences, reverse complements, and
advanced composition calculations with special modes. 
"""
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

    """
    Compute the advanced estimated (probabilistic) composition of the sequence.
    'truncate' mode ignores all special bases, and 'standard' includes all
    special bases without probabilistic conversion.

    Special bases, such as N (which can be any base) or S (which can
    be C or G) for example, can be used to estimate the composition by 
    distributing the occurrence along the possible constituents.
    For example, if the sequence contains 5 S bases, then there's a 
    possibility of a mixture of 5 C or G bases, leading to an estimated
    composition of 2.5 C and 2.5 G.

    This summation uses some advanced short-circuit methods for Python
    total_comp is a dictionary containing the accumulating sum for each
    of the standard bases. The first reduce() iterates through each of
    the special bases in the sequence.

    The second reduce() iterates through each of the possible consituent
    for each special base (i.e. iterate through C and G for S), and adds 
    to the total_comp, generating a new dictionary which becomes pseudo.

    This nested structure performs a reduction along each "axis", one
    being along the special bases, and another along the consituents for
    each special base. Although this is identical to a nested for-loop
    structure, the complications arise with having to accumulate the 
    same data structure for each reduce().

    Parameters
    ----------
    filename : mode (optional)
        Mode of composition computation

    Returns
    -------
    dict[str, int]
        Dictionary of pairs of bases and its occurrence
    """
    def composition (self, mode='estimate') -> dict[str, int]:
        # Generate a Counter that counts all the occurrences of each base
        # len(self.sequence) will return the same as sum(Counter.values())
        occurrence = Counter(self.sequence)

        # Remove all skips
        del occurrence['-']

        # Replace T with U before composition calculation
        # This will be accounted for at the end with uracil presence
        if 'T' in occurrence:
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

        # Estimate the composition
        estimated_comp = dict(reduce(
            lambda total_comp, key_value:
                reduce(
                    lambda pseudo, base:
                        (None if base not in total_comp else pseudo.update({
                            base: total_comp[base] + (
                                key_value[1] / len(SEQ_REPR[key_value[0]])
                                if base in SEQ_REPR[key_value[0]] else 0
                            )
                        })) or pseudo,
                    SEQ_REPR[key_value[0]],
                    total_comp
                ),
            special_bases.items(),
            standard_bases
        ))

        # Replace U with T if DNA
        if not self.uracil_present and 'U' in estimated_comp:
            estimated_comp['T'] = estimated_comp['U']
            del estimated_comp['U']

        return estimated_comp

    def get_name (self) -> str:
        return self.name

    def __repr__ (self) -> str:
        return "<name=\"{}\", sequence={}>".format(self.name, self.sequence)