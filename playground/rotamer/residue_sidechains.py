import networkx as nx
import copy
import matplotlib.pyplot as plt

def add_hydrogens(graph, number_dict):
    """ Adds hydrogens to heavy atoms in a sidechain graph.

    :param graph: sidechain graph
    :param number_dict: dictionary of number of hydrogens per heavy atom
    """
    hydrogens = [''.join(["H", str(i)]) for i in range(0, sum(number_dict.values()))]
    # hydrogens = [atom for atom in graph.nodes() if atom.rstrip("0123456789") == 'H']
    for heavy_atom in number_dict:
        for i in range(0, number_dict[heavy_atom]):
            graph.add_edge(heavy_atom, hydrogens.pop(0))

# A - Alanine
# Not needed

# C - Cysteine
backbone = [("C0", "C1"),
            ("C1", "S0")]
hydrogen_count = {"C1": 2,
                  "S0": 1}
cysteine = nx.Graph(backbone)
add_hydrogens(cysteine, hydrogen_count)

# D - Aspartate
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "O0"),
            ("C2", "O1")]
hydrogen_count = {"C1": 2}
aspartate = nx.Graph(backbone)
add_hydrogens(aspartate, hydrogen_count)

# E - Glutamate
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "C3"),
            ("C3", "O0"),
            ("C3", "O1")]
hydrogen_count = {"C1": 2,
                  "C2": 2}
glutamate = nx.Graph(backbone)
add_hydrogens(glutamate, hydrogen_count)

# F - Phenylalanine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "C3"),
            ("C3", "C4"),
            ("C4", "C5"),
            ("C5", "C6"),
            ("C6", "C7"),
            ("C7", "C2")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "C4": 1,
                  "C5": 1,
                  "C6": 1,
                  "C7": 1}
phenylalanine = nx.Graph(backbone)
add_hydrogens(phenylalanine, hydrogen_count)

# G - Glycine
# Not needed

# Hd - Histidine (delta protonated)

# He - Histidine (epsilon protonated)

# Hp - Histidine (both protonated)

# I - Isoleucine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C1", "C3"),
            ("C3", "C4")]
hydrogen_count = {"C1": 1,
                  "C2": 3,
                  "C3": 2,
                  "C4": 3}
isoleucine = nx.Graph(backbone)
add_hydrogens(isoleucine, hydrogen_count)

# K - Lysine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "C3"),
            ("C3", "C4"),
            ("C4", "N0")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "C3": 2,
                  "C4": 2,
                  "N0": 3}
lysine = nx.Graph(backbone)
add_hydrogens(lysine, hydrogen_count)

# L - Leucine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "C3"),
            ("C2", "C4")]
hydrogen_count = {"C1": 2,
                  "C2": 1,
                  "C3": 3,
                  "C4": 3}
leucine = nx.Graph(backbone)
add_hydrogens(leucine, hydrogen_count)

# M - Methionine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "S0"),
            ("S0", "C3")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "C3": 3}
methionine = nx.Graph(backbone)
add_hydrogens(methionine, hydrogen_count)

# N - Asparagine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "O0"),
            ("C2", "N0")]
hydrogen_count = {"C1": 2,
                  "N0": 2}
asparagine = nx.Graph(backbone)
add_hydrogens(asparagine, hydrogen_count)

# P - Proline
# Not needed

# Q - Glutamine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "C3"),
            ("C3", "O0"),
            ("C3", "N0")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "N0": 2}
glutamine = nx.Graph(backbone)
add_hydrogens(glutamine, hydrogen_count)

# R - Arginine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "C3"),
            ("C3", "C4"),
            ("C4", "N0"),
            ("N0", "C5"),
            ("C5", "N1"),
            ("C5", "N2")]
hydrogen_count = {"C1": 2,
                  "C2": 2,
                  "C3": 2,
                  "C4": 2,
                  "N0": 1,
                  "N1": 2,
                  "N2": 2}
arginine = nx.Graph(backbone)
add_hydrogens(arginine, hydrogen_count)

# S - Serine
backbone = [("C0", "C1"),
            ("C1", "O0")]
hydrogen_count = {"C1": 2,
                  "O0": 1}
serine = nx.Graph(backbone)
add_hydrogens(serine, hydrogen_count)

# T - Threonine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C1", "O0")]
hydrogen_count = {"C1": 1,
                  "C2": 3,
                  "O0": 1}
threonine = nx.Graph(backbone)
add_hydrogens(threonine, hydrogen_count)

# V - Valine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C1", "C3")]
hydrogen_count = {"C1": 1,
                  "C2": 3,
                  "C3": 3}
valine = nx.Graph(backbone)
add_hydrogens(valine, hydrogen_count)

# W - Tryptophan

# Y - Tyrosine
backbone = [("C0", "C1"),
            ("C1", "C2"),
            ("C2", "C3"),
            ("C3", "C4"),
            ("C4", "C5"),
            ("C5", "C6"),
            ("C6", "C7"),
            ("C7", "C2"),
            ("C5", "O0")]
hydrogen_count = {"C1": 2,
                  "C3": 1,
                  "C4": 1,
                  "O0": 1,
                  "C6": 1,
                  "C7": 1}
tyrosine = nx.Graph(backbone)
add_hydrogens(tyrosine, hydrogen_count)

nx.draw_networkx(tyrosine)
plt.show()