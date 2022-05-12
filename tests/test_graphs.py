import sys
sys.path.append('./rosalind')

import unittest

from solutions import graphs

class TestGraphs(unittest.TestCase):
    
    def setUp(self) -> None:
        
        self.overlap_graph_case = './rosalind/problems/overlap_graph_case1.fasta'
        self.answer = ['Rosalind_0498 Rosalind_2391',
                       'Rosalind_0498 Rosalind_0442',
                       'Rosalind_2391 Rosalind_2323']

    def test_find_chain_overlap(self) -> None:
        
        self.assertEqual(self.answer, 
                    graphs.find_chain_overlap(self.overlap_graph_case))

if __name__ == "__main__":
    unittest.main()

