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
        

if __name__ == '__main__':
    def main():
        ## test script
        solver = RosalindSolver()
        # solver.count_nuc('ACGTACGT')
        
        # dna = 'CTTCGAAGTTCATGAGATTCTTCGGACGCCCTCGCTTCATTCGGGCGAGTGTAGAAGAGACGCGCGTGAGATGCTGTAACATAAGATTACGCCTTGTCTGGGTAAGACAGGCGCTATTGGTGGCTGCCGCCGATGGTGTGTACCATCCGTCTTGGAACAATCTGTACAACGTTAACTGTAAGCGGGAATGCCCTAAACGGCGTGACCACTCACATAGACCAGAACGTCTTGCTGGGCCGACCCAGGATCGCTGAGATACGTGAGACGTACGGGAGTAGCACTCTTACGCCGTTGGAAACCATCATTTGCTTGATAACTCGATTGCAATGCGCAGACGCTATTCTCCCACGACCTGATACTAGTGACGTCAGGAGTTAGCTAATAAGCTTGGGTTTGATGTACTTGGGACATGTCTTTAGATTCATTTATCCATGGCTGTCGCATAAAAGTGGTTTTTAACTGAGTAGACTAAGACGTGCGGGGTTAAGAATTTTTCGCCAAAGATGCCGCGAAGCTCGGTGAGACAATCGAGAGTCGTGTGCTGGAACGTAAGCAGCGGCCCACCTGGCGTCATCGACAGGGCGGTCATGCAACCCGTGTTCCTATATACCCGCGTGGCAGACGAGCGTTATCATACTACTTATAACTGTGGATTTACGACTTGCAGAGATTCTGGAGTCGTTCAATACCTATTGATGCGATAAGACAGCTAGATGCAGTGAACCATATTGTACTGGTAAACCTTCTGCACTTAGGATACAGCACCCAAATTACCTTAGACAATCTTACATCTAAATGCATAGCTTTAAATGTATCGCGTGTTTAGATCATGTAGTGGCCAACCCAAGACATCAAGGTAGACACGCTCGCTCGCGCACCTTCGCCCGACACCAACCCTGAGATCTCAGGGCAGGAAGCAAGACGAACGCCTTTTAC'
        # solver.transcribe(dna)
        
        dna = 'ATTCGTGGCTCTGGGGCCCGCGGATAACTGTAATGGCGAAATTGCGGACTATGACCCTGTTATCTAATCACAAAAACGGCGCTAGAAGTGACCCAGAATGTGTGCTGATCCGAATACATCTCACAACAAGTTTACCGCAACGCAACGGGCTTTGCGCTTTTCTAATGATTTGAAGACCGTGGCGAACATTGGCCTAATTACCCACTCTAGTTAATCCCAGACACTGGGGTCTCCAGCGACAGTAAGTCCAGTACATGAGACCAATCTACCAAGTGGTTGGGCACGGCGGTGGAATTACCTCCTCCTTGATCATTACTTTTGACTAATAGCTGTTCGAATTGTAAGCACGCAGGAGTGGCCTGGACGGAGTGTCCGATGCAGTTCGTGACGCCCCTTACGTACTAAGATCGTGGTTAATTCGTGGGATTTCAAAGGAATGGTGTGGTTCTTAGAACTGTCCCAGCCACGTTGGAGGTGCATCTGAATCTCTAGACGGGCCATATTCGGCCAATAGGAATAATCCCCGGGCTTCCTTCTTTAAAGCCACAAGGAGTTAGTAAGGGGGGAGCTATGATCGGTAAGTACGGCGTCCCCGAGGTGGTATCGAAACCGGATGTAACCATATCATCATTGTTATTGACATACGCTGTACTTTATACCAATCTCTCTTGCTCTAAGCGGTGATGTATAGTATATCGCCTCATTCACTGCTGTGACACGGGGAACCTTTGCCAAGCGTACGGTACCTGTAGGAAGTCGAACGGCGTACCGCAGAGCCCAGATCTACGTGCTAGGGCAAGGATATTCCCAGGTTACCAGCAGCATTAAGACTCCGGTATCGCAGAGTCTGCTTTGTCCAACCCGTCCTCGAATTGGTCAAGGCTTGGGTCGCTAGCTCTCGGACGTGTAGCCACTCCCTTTAGACTTTCCAGGTTCCAGCTGCCACGCAAAGCTCACCTCGCCAACCCA'
        solver.complement_strand(dna)

    main()