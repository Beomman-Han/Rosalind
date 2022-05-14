
from typing import List


def calculate_expected_offspring(
    num_pairs : List[int]
    ) -> float:
    
    probs_dominat_phenotype : List[float] = [
        1, 1, 1, .75, .5, 0
    ]
    
    expected_dominant_offspring = 0
    for num, prob in zip(num_pairs, probs_dominat_phenotype):
        expected_dominant_offspring += (2 * num * prob)
    
    return expected_dominant_offspring

if __name__ == '__main__':
    ## solution
    print(calculate_expected_offspring([19785, 16489, 18613, 19160, 17014, 19831]))