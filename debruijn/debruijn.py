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


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return(parser.parse_args())


def read_fastq(fichier) :
    """
    prend un seul argument correspondant au fichier fastq 
    et retourne un générateur de séquences
    """
     
    with open(fichier) as filin:
        for l in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)



def cut_kmer(sequence, k) :
    """
     prend une séquence, une taille de k-mer et retourne un générateur de k-mer
    """
    for i in range(len(sequence)-k+1) :
        yield sequence[i:i+k]


def build_kmer_dict(fichier, k) :
    """
    prend un fichier fastq, une taille k- mer et retourne un dictionnaire 
    ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer
    """
    print("k-mer size : {}".format(k))
    global dico
    dico = {}
    
    for i in read_fastq(fichier) :
        for j in cut_kmer(i, k) :
            if j in dico :
                dico[j] += 1
            else :
                dico[j] = 1           
    print("build_kmer_dict done")
    return(dico)
    

def build_graph(dic) :
    """
    prend en entrée un dictionnaire de k-mer et créera l’arbre 
    de k-mers préfixes et suffixes décrit précédemment
    """
    tree = nx.DiGraph()
    
    for i in dic :
        tree.add_edge(i[:-1], i[1:], weight = dico[i])
    print("build_graph done")
    return tree


def get_starting_nodes(tree) :
    """
    prend en entrée un graphe et retourne une liste de noeuds d’entrée
    """
    starting_node = []
    for i in tree.nodes:
        if len(list(tree.predecessors(i))) == 0 :
            starting_node.append(i)
    print("get_starting_nodes done")
    return starting_node


def get_sink_nodes(tree) :
    """
    prend en entrée un graphe et retourne une liste de noeuds de sortie
    """
    
    sink_nodes = []
    for node in tree.nodes :
        if len(list(tree.successors(node))) == 0 :
            sink_nodes.append(node)
    print("get_sink_nodes done")
    return sink_nodes


def get_contigs(tree, start, sink) :
    """
    prend un graphe, une liste de noeuds d’entrée et une liste de sortie 
    et retourne une liste de tuple(contig, taille du contig)
    """
    
    contigs_list = []

    for i in start :
        for j in sink :
            for l in nx.all_simple_paths(tree, source = i, target = j) :
                c = ""

                for kmer in l :
                    if not c :
                        c += kmer
                    else :
                        c += kmer[-1]
                contigs_list.append((c, len(c)))
    print("get_contigs done")
    return contigs_list


def save_contigs(contig_list, contig_filename) :
    """
    prend une liste de tuple (contig, taille du contig) et un nom de
    fichier de sortie et écrit un fichier de sortie contenant les contigs selon le format fasta
    (retour chariot tous les 80 charactères) à l’aide de la fonction fil
    """
    
    with open(contig_filename, "w") as filout :
        for i, j in enumerate(contig_list) :
            filout.write(f">contig_{i} len={j[1]}\n")
            filout.write(fill(j[0]) + "\n")
    print("save_contigs done")


def fill(text, width = 80) :
    """Split text with a line return to respect fasta format"""
    print("fill done")
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def std(liste : list) :
    """
    prend une liste de valeur, qui retourne l’écart type
    """
    print("std done")
    return statistics.stdev(liste)


def path_average_weight(tree, path) :
    """
     prend un graphe et un chemin et qui retourne un poids
    moyen
    """
    
    weight = 0
    
    for i in range(len(path) - 1) :
        weight += tree.edges[path[i], path[i+1]]["weight"]
    print("path done")
    return weight / (len(path) - 1)


def remove_paths(tree, path_list, delete_entry_node, delete_sink_node) :
    """
    prend un graphe et une liste de chemin, la variable booléenne
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés et la variable
    booléenne delete_sink_node pour indiquer si les noeuds de sortie seront supprimés et
    retourne un graphe nettoyé des chemins indésirables. On détectera une variation entre
    les poids ou les longueurs de nos chemins à l’aide du calcul de l’écart type (std). S’il est
    supérieur à 0, c’est que nous avons des différences
    """
    for path in path_list :
        tree.remove_nodes_from(path[1:-1])
        if delete_entry_node :
            tree.remove_node(path[0])
        if delete_sink_node :
            tree.remove_node(path[-1])
    print("remove_paths done")
    return tree


def select_best_path(tree, path_list, path_length, path_weight, delete_entry_node = False, delete_sink_node = False) :
    """
    prend un graphe, une liste de chemin, une liste donnant la
    longueur de chaque chemin, une liste donnant le poids moyen de chaque chemin,
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés et
    delete_sink_node pour indiquer si les noeuds de sortie seront supprimés et retourne un
    graphe nettoyé des chemins indésirables. Par défaut, delete_entry_node et
    delete_sink_node seront ici à False
    """
    
    path_list_wmax = []
    weight_wmax = []
    length_wmax = []

    for i, j in enumerate(path_list) :
        if path_weight[i] == max(path_weight) :
            path_list_wmax.append(j)
            weight_wmax.append(path_weight[i])
            length_wmax.append(path_length[i])

    path_list_lmax = []
    weight_lmax = []
    length_lmax = []

    for i, j in enumerate(path_list_wmax) :
        if length_wmax[i] == max(length_wmax) :
            path_list_lmax.append(j)
            weight_lmax.append(weight_wmax[i])
            length_lmax.append(length_wmax[i])

    while len(path_list_lmax) > 1 :
        path_list_lmax.pop(statistics.randint(path_list_lmax))

    wrong_paths = []

    for i in path_list :
        if i not in path_list_lmax :
            wrong_paths.append(i)

    tree = remove_paths(tree, wrong_paths, delete_entry_node, delete_sink_node)
    print("select_best_path done")
    return tree


def solve_bubble(tree, a, d) :
    """
    prend un graphe, un noeud ancêtre, un noeud descendant et
    retourne un graph nettoyé de la bulle se trouvant entre ces deux noeuds en utilisant les
    fonctions précédemment développée
    """
    
    liste = []
    path_length = []
    path_weight = []

    for i in nx.all_simple_paths(tree, source = a, target = d) :
        liste.append(i)
        path_length.append(len(i))
        path_weight.append(path_average_weight(tree, i))

    tree = select_best_path(tree, liste, path_length, path_weight, delete_entry_node = False, delete_sink_node = False)
    print("solve_bubble done")
    return tree


def simplify_bubbles(tree) :
    """
    prend un graphe et retourne un graphe sans bulle. Il faut ici
    identifier le noeud ancêtre et le noeud descendant que l’on va 
    indiquer à solve_bubble
    """
    
    wrong_nodes = []

    for i in tree.nodes :
        pred_list = list(tree.predecessors(i))
        if len(pred_list) > 1 :
            a = nx.lowest_common_ancestor(tree, pred_list[0], pred_list[1])
            wrong_nodes.append([a, i])

    for i in wrong_nodes :
        tree = solve_bubble(tree, i[0], i[1])
    print("simplify_bubbles done")
    return tree


def solve_entry_tips(tree, entry) :
    """
    prend un graphe et une liste de noeuds d’entrée et retourne
    graphe sans chemin d’entrée indésirable
    """
    
    ancestors = []
    
    for i in entry :
        for j in nx.descendants(tree, i) :
            if len(tree.pred[j]) >= 2 and j not in ancestors :
                ancestors.append(j)

    paths_list = []
    path_length = []
    path_weight = []

    for i in entry :
        for j in ancestors:
            for l in nx.all_simple_paths(tree, i, j) :
                    paths_list.append(l)
                    path_length.append(len(l))
                    path_weight.append(path_average_weight(tree, l))

        tree = select_best_path(tree, paths_list, path_length, path_weight, delete_entry_node = True, delete_sink_node = False)
        print("solve_entry_tips done")
    return tree


def solve_out_tips(tree, entry) :
    """
    prend un graphe et une liste de noeuds de sortie et retourne
    graphe sans chemin de sortie indésirable
    """
    
    des = []
    
    for i in entry :
        for j in nx.ancestors(tree, i) :
            if len(tree.succ[j]) >= 2 and j not in des :
                des.append(j)

    paths_list = []
    path_length = []
    path_weight = []

    for i in entry :
        for j in des :
            for path in nx.all_simple_paths(tree, j, i) :
                    paths_list.append(path)
                    path_length.append(len(path))
                    path_weight.append(path_average_weight(tree, path))

        tree = select_best_path(tree, paths_list, path_length, path_weight, delete_entry_node = False, delete_sink_node = True)
        print("solve_out_tips done")
    return tree


def main() :
    """
    Main program function
    """
    
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    G = build_graph(kmer_dict)
    start_nodes = get_starting_nodes(G)
    sink_nodes = get_sink_nodes(G)
    contig_list = get_contigs(G, start_nodes, sink_nodes)
    save_contigs(contig_list, args.output_file)
    print("main done")



if __name__ == '__main__':
    main()