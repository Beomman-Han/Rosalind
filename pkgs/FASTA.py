import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from typing import Dict, TextIO, Tuple, Type, Generator

from Seq import Seq
from FileProcessor import FileProcessor
import re, json

class SeqRecord:
    def __init__(
        self,
        seq : str,
        title: str,
        description: str,
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
        self.title = title
        self.description = description
        
        return


class FASTAProcessor(FileProcessor):
    """Class supports functions that process FASTA format file"""
    
    def __init__(self, path : str) -> None:
        """Initialize FASTAProcessor class

        Parameters
        ----------
        path : str
            path of .fasta file
        """
        
        self.path = path
        self.open_obj = False
        
        return
    
    def open(self,
        # file_name : str,
        mode : str = "r"
        ) -> None:
        
        """Open fasta file (self.path) to self.open_obj

        Parameters
        ----------
        # file_name : str
        #     The name of the file to open.
        mode : str
            Open mode (r, w, ...)
        """
        
        if self.open_obj:
            print('Current open_obj is already opened')
        else:
            self.open_mode = mode
            self.open_obj = open(self.path, mode)
        return
    
    def close(self) -> None:
        """Close self.open_obj attribute"""
        self.open_obj.close()
        self.open_obj = False
        return
    
    def readline(self, handle : TextIO) -> Generator[Tuple[str], None, None]:
        """Generator function parsing fasta format contents

        Parameters
        ----------
        handle : TextIO
            TextIO of fasta file

        Yields
        ------
        Tuple[str]
            tuple of contig name, description, sequence
        """
    
        sequences = []
        for line in handle:
            if line.startswith('>'):
                if len(sequences) != 0:
                    yield(title, desc, ''.join(sequences))
                    
                title = line.strip().split()[0][1:]
                try: desc = ' '.join(line.strip().split()[1:])
                except: desc = ''
                sequences = []
            else:
                sequences.append(line.strip())
        yield title, desc, ''.join(sequences)
    
    def write(self,
        title: str,
        sequence: str,
        desc: str = None
        ) -> None:
        
        """Write file by line

        Parameters
        ----------
        title : str
            Title of sequence
        sequence : str
            Sequence
        desc : str
            Description of sequence
        """
        
        if 'r' in self.open_mode:
            print('Current open_obj is "read" mode')        
        elif 'w' in self.open_mode:
            fasta_title = f'>{title}'
            if desc: fasta_title += f' {desc}'
            self.open_obj.write(f'{fasta_title}\n')
            for i in range(0, len(sequence), 70):
                self.open_obj.write(sequence[i:i+70]+'\n')
                
        return
    
    def parse(self, handle : TextIO) -> Generator[Type['SeqRecord'], None, None]:
        """Start parsing fasta file

        Parameters
        ----------
        handle : TextIO
            TextIO of fasta file

        Returns
        -------
        Generator
            self.iterate generator function
        """
        record = self.iterate(handle)
        return record
    
    def iterate(self, handle: TextIO) -> Generator[Type['SeqRecord'], None, None]:
        """Generate function yield SeqRecord object from TextIO

        Parameters
        ----------
        handle : TextIO
            TextIO of fasta file

        Yields
        ------
        SeqRecord
            SeqRecord object containing sequence, title, description
        """
        
        for title, description, seq in self.readline(handle):
            check_out = self.sanity_check((title, description, seq), mode="r")
            if not check_out:
                sys.exit(f'\n==// WARNING //==\n- {title} seqeunce sanity check is failure...\n')
            yield SeqRecord(seq=Seq(seq), title=title, description=description)
            
    def sanity_check(self,
        target_seq : Tuple[str],
        mode : str,
        verbose : bool = False ) -> bool:
        
        """Check if input sequence info ('target_seq') is normal format.

        Parameters
        ----------
        target_seq : Tuple[str]
            Sequence info (title, desc, seq)
        mode : str
            Check mode
        verbose : bool
            Print process message
            (True = Print check message
            False = Silent mode)
            
        Returns
        -------
        bool
            Whether target_seq is normal format
        """

        fasta_elements : list = ['A','T','G','C','N','R','Y','S','W','K','M','B','D','H','V']
        #target_seq = next(fasta_seq)
        if verbose:
            print(f'file format is fasta.\nStart the sanity check for {target_seq[0]}')
    
        if mode == "r":
            check_base = [base.upper() in fasta_elements for base in target_seq[2].strip()]
            check_bool = all(check_base)
            if check_bool == False:
                misbase_list = [ (target_seq[2][i],i+1) for i, b in enumerate(check_base) if b == False]
                #print(misbaseList)
                if verbose:
                    print(f'Please enter a valid seq : \n(misbase,seq position. : {misbase_list} \n')
            if check_bool == False:
                if verbose:
                    print("Check the seq")
                return False
            else:
                if verbose: print("Seq is normal.")
                return True

    def export_to_json(self,
        output_name : str,
        seq_dict : dict = False
        ) -> None:
        
        """Export fasta contents to json format file.
        If 'seq_dict' param is False, then export with self.path fasta info.

        Parameters
        ----------
        output_name : str
            Output file name
        seq_dict : dict, optional
            SeqRecord dictionary, by default False
        """
        
        json_dic = {}
        if not seq_dict:
            for record in self.parse(open(self.path)):
                title = record.title
                json_dic[title] = { 'seq': str(record.seq),
                                    'description': record.description}
        else:
            for title in seq_dict:
                record = seq_dict[title]
                json_dic[title] = { 'seq': str(record.seq),
                                    'description': record.description}
        json.dump(json_dic, open(output_name, 'w'), indent=4)
        return

    def import_from_json(self) -> None:
        pass

    def to_dict(self, handle: TextIO) -> Dict[str, Type['SeqRecord']]:
        """Make dict with {title: SeqRecord object}

        Parameters
        ----------
        handle : TextIO
            TextIO of fasta file

        Returns
        -------
        Dict[str, Type['SeqRecord']]
            dictionary which consists of {title : SeqRecord object}

        Raises
        ------
        ValueError
            Error occurs when duplicate title is in fasta file
        """
        
        record_dict = {}
        for record in self.parse(handle):
            data_id = record.title
            if data_id in record_dict:
                raise ValueError(f"Duplicate key '{data_id}'")
            else:
                record_dict[data_id] = record
        return record_dict

    def find_seq(self, seq: str) -> None:
        fasta_obj = open(self.path,"r")
        p = re.compile(seq)
        for fasta in self.readline(fasta_obj):
            matched_iter = p.finditer((fasta[2]))
            for target in matched_iter:
                print(f'find seq in {fasta[0]} ==> start : {target.start()+1}, end : {target.end()+1}')

    def find_title(self,title:str) -> None:
        fasta_obj = open(self.path,"r")
        p = re.compile(title)
        for fasta in self.readline(fasta_obj):
            matched_iter = p.finditer((fasta[0]))
            for target in matched_iter:
                print(f'find title : {fasta[0]}')
    
if __name__ == "__main__":
    pass