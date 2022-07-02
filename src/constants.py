# https://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
COMPLEMENTARY_BASE = \
{
    # Basic nucleotides for DNA/RNA
    "A": "T",
    "T": "A",
    "U": "A",
    "G": "C",
    "C": "G",
    
    # Special nucleic acid codes from IUPAC
    # Inosine is disregarded as it is not used in the standard IUPAC table,
    # although it pairs with adenine, cytosine, and uracil
    "Y": "R",
    "R": "Y",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",

    # Gap of indeterminable length
    # Since a gap can be any number of unknown nucleotides, its reverse should
    # be a gap of the same length of unknown nucleotides
    "-": "-"
}

# https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#Inverse_RNA_codon_table
# The nucleotide codon uses the compressed form to avoid creating multiple 
# sequences over the same amino acid code
# Arginine, leucine, J-code (leucine/isoleucine), serine, and Z-code (gluatmic
# acid/glutamine) have been compressed as much as possible by hand
INVERSE_AMINO_CODON = \
{
    "A": "GCN",
    "B": "RAY",
    "C": "UGY",
    "D": "GAY",
    "E": "GAR",
    "F": "UUY",
    "G": "GGN",
    "H": "CAY",
    "I": "AUH",
    "J": "NUN",
    "K": "AAR",
    "L": "YUN",
    "M": "AUG",
    "N": "AAY",
    "O": "UAG",
    "P": "CCN",
    "Q": "CAR",
    "R": "MGN",
    "S": "WSN",
    "T": "ACN",
    "U": "UGA",
    "V": "GUN",
    "W": "UGG",
    "X": "NNN",
    "Y": "UAY",
    "Z": "SAR",
    "*": "",
    "-": "-"
}