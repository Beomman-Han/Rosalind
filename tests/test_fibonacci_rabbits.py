import unittest

import sys, os
sys.path.append(os.path.abspath('..'))

from solutions.fibonacci_rabbits import *

class TestFibonacciRabbit(unittest.TestCase):
    """Test ../solutions/fibonacci_rabbits.py"""
    
    def setUp(self) -> None:
        self.input_mortal_rabbits = [(11, 3), (85, 18)]
        self.answer_mortal_rabbits = [16, 258314806822396236]
        
        self.input_immortal_rabbits = [(9, 3)]
        self.answer_immortal_rabbits = [508]
        
    def test_cal_mortal_rabbit(self):
        for input, answer in zip(self.input_mortal_rabbits,
                                self.answer_mortal_rabbits):
            result = cal_mortal_rabbit(*input)
            self.assertEqual(result, answer)
    
    def test_cal_rabbit_pairs(self):
        for input, answer in zip(self.input_immortal_rabbits,
                                self.answer_immortal_rabbits):
            result = cal_rabbit_pairs(*input)
            self.assertEqual(result, answer)
    
    
if __name__ == "__main__":
    unittest.main()