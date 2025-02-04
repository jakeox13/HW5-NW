# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    # Correct Matrices
    m=np.matrix([[0,float("-inf"),float("-inf"),float("-inf"),float("-inf")],
               [float("-inf"), 5., -12., -12., -14.],
               [float("-inf"),-11.,   4.,  -1.,  -6.],
               [float("-inf"),-13.,  -8.,   5.,   4.]])
    ea=np.matrix([[-10., -11., -12., -13., -14.],
                [float("-inf"), -12., -13., -14., -15.],
                [float("-inf"),  -6., -14., -15., -16.],
                [float("-inf"), -7.,  -7., -12., -17.]])
    eb= np.matrix([[-10., float("-inf"), float("-inf"), float("-inf"), float("-inf")],
                [-11., -12.,  -6.,  -7.,  -8.],
                [-12., -13., -14.,  -7.,  -8.],
                [-13., -14., -15., -16.,  -6.]])

    test=NeedlemanWunsch("substitution_matrices/BLOSUM62.mat" , -10,-1)
    test.align(seq1,seq2)
    assert np.array_equal(test._m, m)
    assert np.array_equal(test._ea, ea)
    assert np.array_equal(test._eb, eb)


    pass
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    test=NeedlemanWunsch("substitution_matrices/BLOSUM62.mat" , -10,-1)
    assert test.align(seq3,seq4) == (17,"MAVHQLIRRP","M---QLIRHP")
    pass
