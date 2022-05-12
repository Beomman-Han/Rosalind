"""This module provides class and methods
to deal with sequence data. This mimics
Biopython's Seq and SeqIO object."""

from typing import Generator

#import os, sys, re
#from Seq import Seq


class SeqRecord:
    def __init__(
        self,
        id : str,
        desc : str,
        seq : str
        ) -> None:
        
        """SeqRecord class for recording sequence info (FASTA)

        Parameters
        ----------
        seq : str
            Sequence information
        title : str
            Title of sequence ex) '>title'
        description : str
            Description of sequence
        """
        
        self.seq = seq
        self.id = id
        self.desc = desc
        
        return
    
    @property
    def seq(self) -> str:
        return self._seq.replace('\n', '')
    
    @seq.setter
    def seq(self, _seq : str) -> None:
        self._seq = _seq
    

def parse(fasta : str) -> Generator[SeqRecord, None, None]:
    
    records = []
    for line in open(fasta, 'r'):
        if line.startswith('>'):
            if records:
                yield SeqRecord(*records)
                records = []
        
            labels = line[1:].strip('\n').split()
            ## no description, append empty string
            if len(labels) == 1:
                records += labels
                records += ['']
            else:
                records += [labels[0]]
                ## if description has white space,
                ## bring back to string.
                records += [' '.join(labels[1:])]
        else:
            ## fasta with multi-line sequence
            if len(records) == 3:
                records[-1] += line
            else:
                records.append(line)

    yield SeqRecord(*records)


if __name__ == "__main__":
    pass