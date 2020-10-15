#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
print("Dependencies imported")

print(os.getcwd())

fq = "TCAGAGCTCTAGAGTTGGTTCTGAGAGAGATCGGTTACTCGGAGGAGGCTGTGTCACTCATAGAAGGGATCAATCACACCCACCACGTGTACCGAAACAA"

def read_fastq(fichier):
    global seq
    seq = []
    
    with open(fichier, "r") as filin:
        for line in filin:
            if line.startswith("@"):
                continue
            if line.startswith("+"):
                continue
            if line.startswith("J"):
                continue
            else:
                line = line.replace("\n", "")
                seq.append(line)
    return(seq)


# read_fastq(fq)
# print(seq)

def cut_kmer(sequence, k):
    global kmers
    kmers = []
    l = k
    iii = []
    m = 0

    for i,j in enumerate(sequence):
        if i == k-1:
            for m in range(1,l):
                iii = sequence[i-(l-m)]
            kmers.append(iii)
            k += l
    print(kmers)
    return(kmers)

cut_kmer(fq, 3)

def build_kmer_dict(fichier, k):
    print("k = {}".format(k))
    liste = []
    flattened = []
    read_fastq(fichier)
    
    for i in seq:
        cut_kmer(i, k)
        liste.append([kmers[i] for i in range(len(kmers))])
    
    for sublist in liste:
        for val in sublist:
            flattened.append(val)
    # print(flattened)
    dico = {i:flattened.count(i) for i in flattened}
    print(dico)
    return(dico)
    

# def build_graph(dic):
    


# if __name__ == '__main__':
    # build_kmer_dict(sys.argv[1], int(sys.argv[2]))

