import unittest
from os.path import join, dirname

from src.input_processor import process_fasta_file, process_file_ext, remove_comments



TEST_FILENAME = join(dirname(__file__), "test_file.faa")



class TestInputProcessor (unittest.TestCase):

    def test_file_exists (self):
        with self.assertRaises(FileNotFoundError):
            process_fasta_file('')

    def test_file_ext (self):
        self.assertEqual(process_file_ext(TEST_FILENAME), ".faa")

    def test_comment_removal (self):
        self.assertEqual(
            remove_comments(";comment\n;another\ninformation\n; more comments"),
            "information"
        )

    def test_amino_acid_detection (self):
        sequences = process_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].amino_acid, True)



class TestNucleotideSequence (unittest.TestCase):
    def test_amino_acid_conversion (self):
        sequences = process_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].sequence, "GCN-SAR")

    def test_reverse_complement (self):
        sequences = process_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].reverse_complement().sequence, "YTS-NGC")

if __name__ == "__main__":
    unittest.main()