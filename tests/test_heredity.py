import unittest

import sys
sys.path.append('./rosalind')

from solutions import heredity

class TestDominantOffspring(unittest.TestCase):
    
    def test_calculate_dominant_offspring(self):
        ## Test case 1
        answer = 3.5
        input_ = [1, 0, 0, 1, 0, 1]
        result = heredity.calculate_expected_offspring(input_)
        self.assertAlmostEqual(answer, result)
        
        ## Test case 2
        answer = 13.5
        input_ = [5, 1, 0, 1, 0, 1]
        result = heredity.calculate_expected_offspring(input_)
        self.assertAlmostEqual(answer, result)

    
if __name__ == "__main__":
    unittest.main()