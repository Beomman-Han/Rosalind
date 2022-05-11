"""This is for solving fibonacci rabbits problems.
There are two version of fibonacci rabbits problem,
the one is that rabbits live forever, the other is 
that rabbits only live some of months. Each problem
has each recurrence relation, and both problems
could be solved with 'dynamic programming'.

cal_rabbit_pairs() and cal_rabbit_pairs2() functions
solve the first version (immortal), the former with
tabulation and the latter with memoization. 
cal_mortal_rabbits() function solves the second version (mortal)
with tabulation.
"""

from typing import Dict


def cal_rabbit_pairs(
    n : int,
    k : int
    ) -> int:
    
    """Predict # of rabbit pairs at n th generation in case of reproducing k offsprings
    (implemented by tabulation, bottom-up).

    Recurrence Relation
        Fn = Fn-1 + 3 * Fn-2

    Parameters
    ----------
    n : int
        generation number
    k : int
        offspring number from a reproducable pair

    Returns
    -------
    int
        # of rabbit pairs
    """
    
    ## tabulation
    prev_rabbits, curr_rabbits = 1, 1
    for _ in range(n - 2):
        prev_rabbits, curr_rabbits = curr_rabbits, curr_rabbits + k * prev_rabbits
    
    print(f'[Input]: {n} th generation, {k} offspring from a pair')
    print(f'[Output]: {curr_rabbits} pairs')
    
    return curr_rabbits

def _cal_rabbit_pairs_v2(
    n : int,
    k : int,
    cache : Dict[int, int]
    ) -> int:
    
    """Predict # of rabbit pairs at n th generation in case of reproducing k offsprings
    (implemented by memoization, top-down).

    Recurrence Relation
        Fn = Fn-1 + 3 * Fn-2

    Parameters
    ----------
    n : int
        generation number
    k : int
        offspring number from a reproducable pair
    cache : Dict[int, int]
        Storage containing previous solutions (key: n, value: num_pairs)

    Returns
    -------
    int
        # of rabbit pairs
    """

    ## find from cache
    if n in cache:
        return cache[n]
    
    ## base case
    if n < 3:
        cache[n] = 1
        return cache[n]
    
    ## recursive case
    return k * _cal_rabbit_pairs_v2(n - 2, k, cache) + _cal_rabbit_pairs_v2(n - 1, k, cache)

def cal_rabbit_pairs_v2(
    n : int,
    k : int
    ) -> int:
    
    """Print input, output of _cal_rabbit_pairs_v2.
    
    Parameters
    ----------
    n : int
        generation number
    k : int
        offspring number from a reproducable pair
        
    Returns
    -------
    int
        # of rabbit pairs
    """
    
    num_pairs = _cal_rabbit_pairs_v2(n, k, dict())
    print(f'[Input]: {n} th generation, {k} offspring from a pair')
    print(f'[Output]: {num_pairs} pairs')
    
    return num_pairs

def cal_mortal_rabbit(
    n : int,
    m : int
    ) -> int:
    
    """Calculate the number of rabbits at n th month.
    The rabbits die when they live m month. After 1 month from
    birth, they could reproduce 1 pair of offspring.
    
    Recurrence relation :
    
        F(n) = F(n-1) + F(n-2) (F(1)=F(2)=1, n < m+1)
        F(n) = F(n-1) + ... + F(n-m) (n >= m+1)

    Returns
    -------
    int
        The number of rabbits
    """
    
    if n < 3:
        return 1
    
    generations = [1, 1]
    for k in range(2, n):
        if k < m:
            generations.append(sum(generations[-2:]))
        else:
            if m == 1:
                generations.append(generations[-1])
            else:
                generations.append(sum(generations[k-m:k-1]))
    
    return generations[-1]
