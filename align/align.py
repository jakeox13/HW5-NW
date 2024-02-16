# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        # Match matrix
        self._m = np.empty((len(self._seqB)+1,len(self._seqA)+1))

        # Extension matrices
        self._ea = np.empty((len(self._seqB)+1,len(self._seqA)+1))
        self._eb = np.empty((len(self._seqB)+1,len(self._seqA)+1))

        # Establish trace dictionary
        self.trace=dict()

        # Intialize values first value
        self._m[0,0]=0
        self._ea[0,0]=self.gap_open
        self._eb[0,0]=self.gap_open

        # intilize first row
        for i in range(1,len(seqB)+1):
            self._m[i,0]=float("-inf")
            self._ea[i,0]=float("-inf")
            self._eb[i,0]=self._eb[i-1,0]+self.gap_extend
            self.trace["eb{},0".format(i)]=["b_gap"]
            self.trace["m{},0".format(i)]=["inf"]
            self.trace["ea{},0".format(i)]=["inf"]

        #intilize first column
        for j in range(1,len(seqA)+1):
            self._m[0,j]=float("-inf")
            self._eb[0,j]=float("-inf")
            self._ea[0,j]=self._ea[0,j-1]+self.gap_extend
            self.trace["ea0,{}".format(j)]=["a_gap"]
            self.trace["m0,{}".format(j)]=["inf"]
            self.trace["eb0,{}".format(j)]=["inf"]

        self.trace["m0,0"]=["end"]

        # Loop through and fill out matrix
        for i in range(1,len(seqB)+1):
            for j in range(1,len(seqA)+1):
                # Calculate match options
                m_opt=[    self._m[i-1,j-1]+self.sub_dict[(seqB[i-1],seqA[j-1])],
                           self._ea[i-1,j-1]+self.sub_dict[(seqB[i-1],seqA[j-1])],
                           self._eb[i-1,j-1]+self.sub_dict[(seqB[i-1],seqA[j-1])]
                           ]
                           
                self._m[i,j]=max(m_opt)

               
                #Get the index of the max value since it takes first value this priottizes from m
                prev=m_opt.index(max(m_opt))
                if prev==0:
                    self.trace["m{0},{1}".format(i,j)]= self.trace["m{0},{1}".format(i-1,j-1)]+["diag"]
                elif prev==1: 
                    self.trace["m{0},{1}".format(i,j)]= self.trace["ea{0},{1}".format(i-1,j-1)]+["diag"]
                else:
                    self.trace["m{0},{1}".format(i,j)]= self.trace["eb{0},{1}".format(i-1,j-1)]+["diag"]
                
                # Calculate extend a options
                ea_opt=[self._m[i-1,j]+self.gap_open+self.gap_extend,
                        self._ea[i-1,j]+self.gap_extend,
                        self._eb[i-1,j]+self.gap_open+self.gap_extend]

                self._ea[i,j]=max(ea_opt)

                prev=ea_opt.index(max(ea_opt))
                if prev==0:
                    self.trace["ea{0},{1}".format(i,j)]= self.trace["m{0},{1}".format(i-1,j)]+["a_gap"]
                elif prev== 1:
                    self.trace["ea{0},{1}".format(i,j)]= self.trace["ea{0},{1}".format(i-1,j)]+["a_gap"]
                else: 
                    self.trace["ea{0},{1}".format(i,j)]= self.trace["eb{0},{1}".format(i-1,j)]+["a_gap"]
               
                # Calculate extend b options
                eb_opt=[self._m[i,j-1]+self.gap_open+self.gap_extend,
                        self._ea[i,j-1]+self.gap_open+self.gap_extend,
                        self._eb[i,j-1]+self.gap_extend]
                self._eb[i,j]=max(eb_opt)

                prev=eb_opt.index(max(eb_opt))
                if prev==0:
                    self.trace["eb{0},{1}".format(i,j)]= self.trace["m{0},{1}".format(i,j-1)]+["b_gap"]
                elif prev==1:
                    self.trace["eb{0},{1}".format(i,j)]= self.trace["ea{0},{1}".format(i,j-1)]+["b_gap"]
                else: 
                    self.trace["eb{0},{1}".format(i,j)]= self.trace["eb{0},{1}".format(i,j-1)]+["b_gap"]

                # print(m)
                # print(ea)
                # print(eb)

        # pull final values for traceback
        self._m_final = self._m[(len(self._seqB),len(self._seqA))]
        self._ea_final = self._ea[(len(self._seqB),len(self._seqA))]
        self._eb_final = self._eb[(len(self._seqB),len(self._seqA))]	


        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        end_res=[self._m_final,self._ea_final,self._eb_final]
        self.alignment_score=max(end_res)

        # Figure out path
        prev=end_res.index(max(end_res))
        if prev==0:
            path=self.trace["m{},{}".format(len(self._seqB),len(self._seqA))]
        elif prev==1:
            path=self.trace["ea{},{}".format(len(self._seqB),len(self._seqA))]
        else:
            path=self.trace["eb{},{}".format(len(self._seqB),len(self._seqA))]
        
        # Start at front of each sequence
        a_index=0
        b_index=0
        self.seqA_align=""
        self.seqB_align=""

        # Start at 1 since first should always be end
        for val in path[1:]:
            # If this is diagonal 
            if val == "diag":
                # Add the next value from both sequences
                self.seqA_align="".join([self.seqA_align, self._seqA[a_index]])
                self.seqB_align="".join([self.seqB_align, self._seqB[b_index]])
                
                # Move to next potential value
                a_index+=1
                b_index+=1

            elif val == "a_gap":
                self.seqA_align="".join([self.seqA_align, "-"])
                self.seqB_align="".join([self.seqB_align, self._seqB[b_index]])

                b_index+=1

            # b_gap senario
            else:
                self.seqA_align="".join([self.seqA_align, self._seqA[a_index]])
                self.seqB_align="".join([self.seqB_align, "-"])

                a_index+=1




        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
