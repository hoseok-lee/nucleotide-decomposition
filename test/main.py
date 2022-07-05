import unittest
from os.path import join, dirname

import src.fasta_processor as fasta_proc
from src.nucleotide_sequence import NucleotideSequence



TEST_FILENAME = join(dirname(__file__), "test_file.faa")



class TestInputProcessor (unittest.TestCase):
    def test_file_ext (self):
        self.assertEqual(fasta_proc.get_file_ext(TEST_FILENAME), ".faa")

    def test_fasta_chunk_decode_encode (self):
        test_fasta_chunk = {"name": "sequence", "another name": "next sequence"}
        encoded = fasta_proc.encode_fasta_chunks(test_fasta_chunk)
        self.assertEqual(
            fasta_proc.decode_fasta_chunks(encoded),
            test_fasta_chunk
        )

    def test_amino_acid_detection (self):
        sequences = fasta_proc.open_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].amino_acid, True)



class TestNucleotideSequence (unittest.TestCase):
    def test_amino_acid_conversion (self):
        sequences = fasta_proc.open_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].sequence, "GCN-SAR")

    def test_reverse_complement (self):
        sequences = fasta_proc.open_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].reverse_complement().sequence, "YTS-NGC")

    def test_composition (self):
        sequence1 = NucleotideSequence(name="", sequence="ATGCATGC")
        sequence2 = NucleotideSequence(name="", sequence="RYSWNNNN")
        self.assertEqual(sequence1.composition(), sequence2.composition())

    def test_composition_verbose (self):
        sequence1 = NucleotideSequence(name="", sequence="RYNNNNSW")
        sequence2 = NucleotideSequence(name="", sequence="NNRYSWNN")
        self.assertEqual(sequence1.composition(), sequence2.composition())

if __name__ == "__main__":
    unittest.main()