# -*- coding: utf-8 -*-
"""
Sequence Aligment
	- Smith-Waterman algorithm

Xavier Pinho & Jorge Melo- Introduction to Bioinformatics, University of Coimbra, 2018/2019
"""

import numpy as np
import os

print("--------Importation Complete----------")

home = os.getcwd()
os.chdir(home)

file_dna = open(home+"/sequence_dna.txt")
sequence_dna = file_dna.read()
sequence_dna = sequence_dna.replace("\n","")

file_mrna = open(home+"/sequence_mrna.txt")
sequence_mrna = file_mrna.read()
sequence_mrna = sequence_mrna.replace("\n","")

print("-----------Data Loaded------------------")

def smith_waterman(a: str, b: str, alignment_score: float = 1, gap_cost: float = 2) -> float:
    """
    Compute the Smith-Waterman alignment score for two strings.
    See https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Algorithm
    This implementation has a fixed gap cost (i.e. extending a gap is considered
    free). In the terminology of the Wikipedia description, W_k = {c, c, c, ...}.
    This implementation also has a fixed alignment score, awarded if the relevant
    characters are equal.
    Kinda slow, especially for large (50+ char) inputs.
    """
    # H holds the alignment score at each point, computed incrementally
    H = np.zeros((len(a) + 1, len(b) + 1))
    P = np.zeros((len(a) + 1, len(b) + 1))
    DELETION, INSERTION, MISMATCH, MATCH = range(4)

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            # The score for substituting the letter a[i-1] for b[j-1]. Generally low
            # for mismatch, high for match.
            match    = (H[i - 1, j - 1] + (alignment_score if a[i - 1] == b[j - 1] else 0), MATCH)
            mismatch = (H[i - 1, j - 1] - (alignment_score if a[i - 1] != b[j - 1] else 0), MISMATCH)
            # The scores for for introducing extra letters in one of the strings (or
            # by symmetry, deleting them from the other).
            delete = (H[1:i, j].max() - gap_cost if i > 1 else 0, DELETION)
            insert = (H[i, 1:j].max() - gap_cost if j > 1 else 0, INSERTION)
            H[i, j], P[i, j] = max(match, delete, insert, mismatch, (0, 0))

    outfile1 = "part1_database"
    np.save(outfile1,H)
    outfile2 = "part1_database_p"
    np.save(outfile2,P)

smith_waterman(a=sequence_dna, b=sequence_mrna)
print("---------Done---------------")
