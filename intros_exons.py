# -*- coding: utf-8 -*-
"""
Sequence Aligment
    - Finding exons and introns

Xavier Pinho & Jorge Melo- Introduction to Bioinformatics, University of Coimbra, 2018/2019
"""
import numpy as np
import os

print("--------Importation Complete----------")

home = os.getcwd()
os.chdir(home)

P = np.load(home+"/database_p.npy")
H = np.load(home+"/database.npy")

file_dna = open(home+"/sequence_dna.txt")
sequence_dna = file_dna.read()
sequence_dna = sequence_dna.replace("\n","")

file_mrna = open(home+"/sequence_mrna.txt")
sequence_mrna = file_mrna.read()
sequence_mrna = sequence_mrna.replace("\n","")

print("-----------Data Loaded------------------")

ind = np.unravel_index(np.argmax(H, axis=None), H.shape)
i, j = ind

DELETION, INSERTION, MISMATCH, MATCH = range(4)

def run():

    def backtrack():
        i, j = ind
        while i > 0 or j > 0:
            assert i >= 0 and j >= 0
            if P[i][j] == MATCH or P[i][j] == MISMATCH:
                i -= 1
                j -= 1
                yield sequence_dna[i], sequence_mrna[j]
            elif P[i][j] == INSERTION:
                j -= 1
                yield '-', sequence_mrna[j]
            elif P[i][j] == DELETION:
                i -= 1
                yield sequence_dna[i], '-'
            else:
                assert (False)


    return [''.join(reversed(s)) for s in zip(*backtrack())]


a = run()
result_dna  = a[0]
result_mrna = a[1]
result_mrna = result_mrna.replace("-","")

exons = [result_dna[0:91], result_dna[223:940], result_dna[1071:1118], result_dna[1219:2260], result_dna[2409:3399]]
introns = [result_dna[92:222], result_dna[941:1070], result_dna[1119:1218], result_dna[2261:2408]]
print("Exons: ", exons)
print("Introns: ",introns)

print("------------Done-------------")