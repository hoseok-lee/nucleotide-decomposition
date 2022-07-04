import unittest
from os.path import join, dirname

import src.input_processor as inp_proc



TEST_FILENAME = join(dirname(__file__), "test_file.faa")



class TestInputProcessor (unittest.TestCase):

    def test_file_exists (self):
        with self.assertRaises(FileNotFoundError):
            inp_proc.process_fasta_file('')

    def test_file_ext (self):
        self.assertEqual(inp_proc.process_file_ext(TEST_FILENAME), ".faa")

    def test_comment_removal (self):
        test_string = ";comment\ninformation\n; more comments"
        self.assertEqual(inp_proc.remove_comments(test_string), "information")

    def test_fasta_chunk_generation (self):
        test_string = ">name\nsequence>another name\nnext sequence"
        self.assertEqual(
            inp_proc.generate_fasta_chunks(test_string),
            {"name": "sequence", "another name": "next sequence"}
        )

    def test_amino_acid_detection (self):
        sequences = inp_proc.process_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].amino_acid, True)



class TestNucleotideSequence (unittest.TestCase):
    def test_amino_acid_conversion (self):
        sequences = inp_proc.process_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].sequence, "GCN-SAR")

    def test_reverse_complement (self):
        sequences = inp_proc.process_fasta_file(TEST_FILENAME)
        self.assertEqual(sequences[0].reverse_complement().sequence, "YTS-NGC")

if __name__ == "__main__":
    unittest.main()