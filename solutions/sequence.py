"""This contains solutions for problems relative to
DNA/RNA/Protein sequence."""

from typing import Dict, List
#from Bio import SeqIO

def find_common_substring(
    fasta : str
    ) -> str:
    
    ## load seequences from fasta file
    sequences : List[str] = []
    fin_stream = open(fasta, 'r')
    sequence = f''
    for line in fin_stream:
        if line.startswith('>'):
            if sequence:
                sequences.append(sequence)
                sequence = f''
            continue
        sequence += line.strip()
    if sequence:
        sequences.append(sequence)
    fin_stream.close()
    
    print(sequences)
    
    if len(sequences) < 2:
        return NotImplemented

    ## firstly, find initial common substrings from 2 sequences
    common_substrings : Dict[int, List[str]] = dict()
    for i, s1 in enumerate(sequences[0]):
        for j, s2 in enumerate(sequences[1]):
            if s1 == s2:
                substring = find_longest_substring(
                    sequences[0][i:], sequences[1][j:])
                if len(substring) > 1:
                    try:
                        common_substrings[len(substring)].append(substring)
                    except KeyError:
                        common_substrings[len(substring)] = [substring]
    
    print(common_substrings)
    
    ## find initial common substrings at other sequences
    # for length in sorted(list(common_substrings.keys()), reverse=True):
    for seq in sequences[2:]:
        new_common_substrings : Dict[int, List[str]] = dict()
        for length in common_substrings.keys():
            if len(seq) >= length:
                for substr in common_substrings[length]:
                    if seq.find(substr) != -1:
                        try:
                            new_common_substrings[length].append(substr)
                        except KeyError:
                            new_common_substrings[length] = [substr]
        common_substrings = new_common_substrings
        print(common_substrings)
    
    if common_substrings:
        return common_substrings[max(common_substrings.keys())][0]
    
    return ''
        
def find_longest_substring(
    seq1 : str,
    seq2 : str
    ) -> str:
    
    longest_substring = f''
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] == seq2[i]:
            longest_substring += f'{seq1[i]}'
        else:
            break

    return longest_substring


if __name__ == '__main__':
    # common_substr : str = find_common_substring('../problems/find_common_substrings_case1.fasta')
    common_substr : str = find_common_substring('../problems/longest_common_substring_case.fasta')
    print(common_substr)