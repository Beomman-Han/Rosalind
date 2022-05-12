"""This is for solving problems using graph data structure.
The basic graph problems could be extended to assembly
algorithm."""

from typing import List

import sys
sys.path.append('./rosalind')

from pkgs import FASTA


def find_chain_overlap(fasta : str) -> List[str]:
    """Find chain overlapped string and make 
    directed graph format, which is 'overlap graph.
    
    For a collection of strings and a positive integer k,
    the overlap graph for the strings is a directed 
    graph O(k) where each string is represented by
    a node, and string s is connected to string t with
    a directed edge if a length k suffix of s matches
    a length k prefix of t (but s != t).
    
    A directed graph contains directed edges which consist
    of tail and head. Thus a directed graph could be
    represented by list of directed edges of tuple (tail, head).
    
    Parameters
    ----------
    fasta : str
        Path of a fasta file which has labeling and sequence

    Returns
    -------
    List[str]
        A directed graph implemented with list of string
    """
    
    overlap_graph : List[str] = []
    for record in FASTA.parse(fasta):
        for record2 in FASTA.parse(fasta):
            if record.id == record2.id:
                continue
            if record.seq[-3:] == record2.seq[:3]:
                overlap_graph.append([record.id, record2.id])
    
    return [ f'{edge[0]} {edge[1]}' for edge in overlap_graph ]

if __name__ == "__main__":
    ## solution
    for edge in find_chain_overlap('./rosalind/problems/overlap_graph_case.fasta'):
        print(edge)
    
