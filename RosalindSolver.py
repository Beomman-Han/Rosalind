from typing import Dict, Tuple


class RosalindSolver:
    """Class for solving 'Rosalind' bioinformatics problems.
    Each method represents each Rosalind problem.
    """
    
    DNA_BASES = ('A', 'C', 'G', 'T')
    DNA_WC_PAIR = {'A': 'T', 'T': 'A',
                   'C': 'G', 'G': 'C'}
    
    def __init__(self):
        pass
    
    def count_nuc(self,
        seq : str
        ) -> Dict[str, int]:
        
        """Count A,C,G,T nucleotide from input DNA sequence.

        Parameters
        ----------
        seq : str
            Input DNA sequence

        Returns
        -------
        Dict[str, int]
            Count of A, C, G, T
        """

        ## init dictionary for counting        
        nuc_counts = {nuc : 0 for nuc in RosalindSolver.DNA_BASES}
        for nuc in seq.upper():
            try:
                nuc_counts[nuc] += 1
            except KeyError:
                continue
        
        print(f'[Input]: {seq}')
        print(f'[Output]:', end=' ')
        for nuc in RosalindSolver.DNA_BASES:
            print(nuc_counts[nuc], end=' ')
        print()

        return nuc_counts
    
    def transcribe(self, dna_seq : str) -> str:
        """Transcribe DNA sequence to RNA sequence.
        (It does not consider the direction of molecule, only replace 'T' to 'U')

        Parameters
        ----------
        dna_seq : str
            Input DNA sequence

        Returns
        -------
        str
            RNA sequence from input DNA
        """
        
        rna_seq = dna_seq.replace('T', 'U').replace('t', 'u')
        print(f'[Input]: {dna_seq}')
        print(f'[Output]: {rna_seq}')
        return rna_seq
    
    def complement_strand(self, seq : str) -> str:
        """Get reverse complementary sequence of input sequence.

        Parameters
        ----------
        seq : str
            Input DNA sequence

        Returns
        -------
        str
            Sequence of complement strand
        """
        
        comp_seq = ''
        for nuc in seq:
            comp_seq += RosalindSolver.DNA_WC_PAIR[nuc]
        comp_seq = comp_seq[::-1]
        
        print(f'[Input]: {seq}')
        print(f'[Output]: {comp_seq}')
        
        return comp_seq
    
    def cal_rabbit_pairs(self,
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
    
    def _cal_rabbit_pairs_v2(self,
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
        return k * self._cal_rabbit_pairs_v2(n - 2, k, cache) + self._cal_rabbit_pairs_v2(n - 1, k, cache)
    
    def cal_rabbit_pairs_v2(self,
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
        
        num_pairs = self._cal_rabbit_pairs_v2(n, k, dict())
        print(f'[Input]: {n} th generation, {k} offspring from a pair')
        print(f'[Output]: {num_pairs} pairs')
        
        return num_pairs

    def get_highest_gc(self,
        fasta : str
        ) -> str:
        
        """Get sequence ID and GC ratio which has the highesst GC(%)
        from input fasta file.

        Parameters
        ----------
        fasta : str
            Path of input fasta file

        Returns
        -------
        str
            Seq ID and GC ratio
        """
        
        ## initialize variables
        id_max_gc, max_gc = '', 0
        
        fasta = open(fasta, 'r')
        ## initialize variables
        seq_id, seq = '', ''
        for line in fasta:
            if line[0] == '>':
                ## sequence id line
                if len(seq) == 0:
                    ## if no seq, reset directly
                    seq_id = line[1:].strip()
                    seq = ''
                    continue
                
                ## save previous sequence information
                gc_count = seq.count('G') + seq.count('C')
                gc_ratio = round(100 * gc_count / len(seq), 6)
                
                if gc_ratio >= max_gc:
                    id_max_gc = seq_id
                    max_gc = gc_ratio
                
                ## reset
                seq_id = line[1:].strip()
                seq = ''
            else:
                seq += line.strip().upper()
        fasta.close()
        
        ## sequence id line
        if len(seq) != 0:
            ## save previous sequence information
            gc_count = seq.count('G') + seq.count('C')
            gc_ratio = round(100 * gc_count / len(seq), 6)
        
            if gc_ratio >= max_gc:
                id_max_gc = seq_id
                max_gc = gc_ratio
        
        answer = f'{id_max_gc}\n{max_gc}'
        print(answer)
        
        return answer
    
    def count_point_mutation(self,
        seq1 : str,
        seq2 : str
        ) -> int:
        
        """Count point mutation between s, t sequence.
        (equals to calculate 'Hamming' distance)

        Parameters
        ----------
        seq1 : str
            Input DNA string 1
        seq2 : str
            Input DNA string 2

        Returns
        -------
        int
            Number of diff nucleotides
        """
        
        hamming_dist = 0
        for s, t in zip(seq1, seq2):
            if s != t:
                hamming_dist += 1
        
        print(hamming_dist)
        
        return hamming_dist
    

if __name__ == '__main__':
    def main():
        ## test script
        solver = RosalindSolver()
        # solver.count_nuc('ACGTACGT')
        
        # dna = 'CTTCGAAGTTCATGAGATTCTTCGGACGCCCTCGCTTCATTCGGGCGAGTGTAGAAGAGACGCGCGTGAGATGCTGTAACATAAGATTACGCCTTGTCTGGGTAAGACAGGCGCTATTGGTGGCTGCCGCCGATGGTGTGTACCATCCGTCTTGGAACAATCTGTACAACGTTAACTGTAAGCGGGAATGCCCTAAACGGCGTGACCACTCACATAGACCAGAACGTCTTGCTGGGCCGACCCAGGATCGCTGAGATACGTGAGACGTACGGGAGTAGCACTCTTACGCCGTTGGAAACCATCATTTGCTTGATAACTCGATTGCAATGCGCAGACGCTATTCTCCCACGACCTGATACTAGTGACGTCAGGAGTTAGCTAATAAGCTTGGGTTTGATGTACTTGGGACATGTCTTTAGATTCATTTATCCATGGCTGTCGCATAAAAGTGGTTTTTAACTGAGTAGACTAAGACGTGCGGGGTTAAGAATTTTTCGCCAAAGATGCCGCGAAGCTCGGTGAGACAATCGAGAGTCGTGTGCTGGAACGTAAGCAGCGGCCCACCTGGCGTCATCGACAGGGCGGTCATGCAACCCGTGTTCCTATATACCCGCGTGGCAGACGAGCGTTATCATACTACTTATAACTGTGGATTTACGACTTGCAGAGATTCTGGAGTCGTTCAATACCTATTGATGCGATAAGACAGCTAGATGCAGTGAACCATATTGTACTGGTAAACCTTCTGCACTTAGGATACAGCACCCAAATTACCTTAGACAATCTTACATCTAAATGCATAGCTTTAAATGTATCGCGTGTTTAGATCATGTAGTGGCCAACCCAAGACATCAAGGTAGACACGCTCGCTCGCGCACCTTCGCCCGACACCAACCCTGAGATCTCAGGGCAGGAAGCAAGACGAACGCCTTTTAC'
        # solver.transcribe(dna)
        
        # dna = 'ATTCGTGGCTCTGGGGCCCGCGGATAACTGTAATGGCGAAATTGCGGACTATGACCCTGTTATCTAATCACAAAAACGGCGCTAGAAGTGACCCAGAATGTGTGCTGATCCGAATACATCTCACAACAAGTTTACCGCAACGCAACGGGCTTTGCGCTTTTCTAATGATTTGAAGACCGTGGCGAACATTGGCCTAATTACCCACTCTAGTTAATCCCAGACACTGGGGTCTCCAGCGACAGTAAGTCCAGTACATGAGACCAATCTACCAAGTGGTTGGGCACGGCGGTGGAATTACCTCCTCCTTGATCATTACTTTTGACTAATAGCTGTTCGAATTGTAAGCACGCAGGAGTGGCCTGGACGGAGTGTCCGATGCAGTTCGTGACGCCCCTTACGTACTAAGATCGTGGTTAATTCGTGGGATTTCAAAGGAATGGTGTGGTTCTTAGAACTGTCCCAGCCACGTTGGAGGTGCATCTGAATCTCTAGACGGGCCATATTCGGCCAATAGGAATAATCCCCGGGCTTCCTTCTTTAAAGCCACAAGGAGTTAGTAAGGGGGGAGCTATGATCGGTAAGTACGGCGTCCCCGAGGTGGTATCGAAACCGGATGTAACCATATCATCATTGTTATTGACATACGCTGTACTTTATACCAATCTCTCTTGCTCTAAGCGGTGATGTATAGTATATCGCCTCATTCACTGCTGTGACACGGGGAACCTTTGCCAAGCGTACGGTACCTGTAGGAAGTCGAACGGCGTACCGCAGAGCCCAGATCTACGTGCTAGGGCAAGGATATTCCCAGGTTACCAGCAGCATTAAGACTCCGGTATCGCAGAGTCTGCTTTGTCCAACCCGTCCTCGAATTGGTCAAGGCTTGGGTCGCTAGCTCTCGGACGTGTAGCCACTCCCTTTAGACTTTCCAGGTTCCAGCTGCCACGCAAAGCTCACCTCGCCAACCCA'
        # solver.complement_strand(dna)
        
        # n, k = 33, 5
        # solver.cal_rabbit_pairs(n, k)
        # solver.cal_rabbit_pairs_v2(n, k)
        
        # solver.get_highest_gc('/Users/hanbeomman/Downloads/rosalind_gc.txt')
        
        seq1 = 'AGGGGAGGCCCGACAAATGCCAAACTTATGCATATGGTTTTTCCGAGCACTACGATCTCCGCTCCAGTGCGAATCTCTGAGTACGAGTATTCCAGGTCAGACTCGCATCCTTAAACGCTCCGATGACCGATTATGTCTAATGTTTAAAGGCAATCGCTGTAGATAACCACCTTTAAAAACCAGACATGGGTACTGCGCAGTCCGTAGAATACCGTTCGTGATAGTGCATGTTTCGGTCCACCTACGGTCATCCGCGTCCGGAGGCAAGGTTATGACCACACCCTGTATCGTAACAGGGTGAATGCCAATAGTAGCGGCACTTCTCATGCCACCCGTGTCCGGCGCTTTGCTATGTACTTCACCGCTTAATTCCAAACGGCAGAGGGGTTCGCCCTGGGAAGGTTCTTCCAATACATCACATTTTGCCAGAGAGTTGGCACTCGTAGGTTGTGTCAAGAGTCCTTTGCGGGTTCGGAGCCGATAAAAGTCATGGCTCCTCTGTATCTTCTCGATCTGGTTCGTTATGAGTATATGATGTGAGTTTGGCATTAAAGTAATTGGTTAGGATTACCTGAAGCGTCCAGATCTTGCATTTCCCACTTACATGAATAATGTCCAATTTTACGGCGCGTACCCATGCATTTACAAGCACGTTTTGCCTTAATTGATATTTCAAAGAAGCCCTAGATGTCACCTAGTAAGGACCTAGCCTTATACACGTTATAAATGTGTTTCTTTGCTTATAAGTACGTCCATCAAATCCTTTCCCTGTGATGCCATCGTCATGACGCTTGATCCGCCCTGTCCTTCGTACATGTACGATCCCGACGTGGCCCTCTTCTTAGGGTGTGCTGCAAATATGGTTCCCGACTGCCCTCCGCGACGAACGAAACGTTTCCCCTCCTAACACTCCAGGTCCCGTGTAAGAACAGCGCGCCCTGAAGCATAATTTGCCACGGCCTTTCGCTTCCCTGTTGCTTATA'
        seq2 = 'AGTGGAGGCGCGACGTGGATCAATGATTTGGATTGGATTGTTAGGTCCTCATGCGAGTGCGCGGGTTGAAGATGCACTAAGGTCCGCTATGCCAGATCAGACTCGGCTTCTTAAACGGTCCGATGACCCCTTAAAGCTTCAGATTAAGGGATAGTGTTTTATATATTCACCTGTACCACCCAATGATGGGTCATGCCGATTCCCTACCTTAACGTTTTATAGAGAGGCCGCTAGGAAGTTACGCCGACCTTCCCAGCTATGAGAAATGATCCAGTAGCGGCCGGTTATCCTCCACGGGTTCTCGGCGGTGGTGTAGGGTTATCAGAAGTCCCTAAAAGCCTTTCGATTGTTGTGTATCTCACTACGCCATCCTATACCTTTACAGACTTGGCGACAGTAACTATCGCCCCGTACATCACTTTGAACAAGGGAGACGGCGATGGTAATCTATGCACATTGTCCCGATCCGGTTCCGTTCCTATGGAGGTCCTGGCGAGTCATTAGATTCGTTATCTGGTTTGTCATGGGTAGACGATTGTACCCGGCTTTTTAATGTTTCCCTATGGTGTACCGCAGGCTGCCAGTTCAGTCGTCTGGTCATTCAATGCATAGAGGGCAAGCTGACCGCGCGCATCTATGTGCTTCATAACCCCTCGGCGCTTCAGTCCTGGTACGGAGAAGACCTAAGATTCACCTACAAGGTACCCAGCCTGAAAATGGTGTTAAGAGAGTATATCAGCTTATAAGTTTGTTCATTAGTTCGTTTCCCTATGCTTTCGTGGTGCTGACGAAAAGTCCGGCCCCTCACCCGAACTAGGATGTGCCCGTCGTGGTCACCTGCATATGGGTAGTAATGCCCGGGTGTTCCGATTAAACGCCGCGATCAACGGTACTGACCCTGTTCCAAAACTGCGATACCACTGACCAATGAGCACACCAGGTAGCAAACTTACCCACGACGGCTCGTTCATGCATTACCTACT'
        solver.count_point_mutation(seq1, seq2)
        
    main()