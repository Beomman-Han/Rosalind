
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

def independent_alleles(
    generation : int,
    least_num : int
    ) -> float:
    
    num_population = 2 ** generation
    prob = 0.0
    for i in range(least_num, num_population + 1):
        prob += _get_factorial(num_population) / \
            (_get_factorial(num_population - i) * _get_factorial(i)) * \
                (0.25 ** i) * (0.75 ** (num_population - i))
    return prob

def _get_factorial(n : int) -> int:
    factorial = 1
    for i in range(1, n + 1):
        factorial *= i
    return factorial
    

if __name__ == '__main__':
    ## solution
    print(calculate_expected_offspring([19785, 16489, 18613, 19160, 17014, 19831]))
    print(independent_alleles(5, 7))