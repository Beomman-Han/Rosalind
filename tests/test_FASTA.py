import os, sys
sys.path.append('./rosalind')

from pkgs import FASTA

import unittest

class TestFASTA(unittest.TestCase):
    
    def setUp(self) -> None:
        self.fasta = './rosalind/problems/overlap_graph_case1.fasta'
        self.ids = ['Rosalind_0498', 'Rosalind_2391',
                    'Rosalind_2323', 'Rosalind_0442']
        self.seqs = ['AAATAAA', 'AAATTTT',
                     'TTTTCCC', 'AAATCCC']
    
    def test_parse(self) -> None:
        for idx, seq_record in enumerate(FASTA.parse(self.fasta)):
            self.assertEqual(seq_record.id, self.ids[idx])
            self.assertEqual(seq_record.seq, self.seqs[idx])
        
if __name__ == "__main__":
    unittest.main()