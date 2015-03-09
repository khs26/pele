#!/usr/bin/python

import my_read_amber as ra
import networkx as nx
import sys

class GhostAtom(ra.Atom):
  pass

valences = { "H": 1,
         "C": 4,
         "N": 3,
         "O": 2,
         "S": 2,
         "Se": 2,
         "P": 5 }

def chiral_candidates(atoms):
  candidates = [atom for atom in atoms.nodes() if len(nx.edges(atoms, atom)) == 4]
  return candidates

def multi_bonds(atoms):
  multibonded = [atom for atom in atoms.nodes() 
                 if len(nx.edges(atoms, atom)) < valences[atom.element]]
  for i, atom in enumerate(multibonded):
    paired = False
    for other in atoms.neighbors(atom):
      if isinstance(other, GhostAtom):
        paired = True
        continue
      if len(nx.edges(atoms, other)) < valences[other.element]:
        ghost_atom = GhostAtom(**(atom.__dict__))
        ghost_atom.name = atom.name + "*"
        ghost_other = GhostAtom(**(other.__dict__))
        ghost_other.name = other.name + "*"
        atoms.add_edge(other, ghost_atom)
        atoms.add_edge(atom, ghost_other)
        paired = True
#    if not paired:
#      print "Could not find pair for atom", atom, atom.residue, atom.charge, nx.edges(atoms, atom)

def chiral_order(atoms, chiral_atom, depth=6):
    # print "\n\nResidue:", chiral_atom.residue, "atom:", chiral_atom
    # print "Neighbours:", atoms.neighbors(chiral_atom)
    # Create a list of ordered atoms to be passed back
    ordered = []
    # Do a quick check whether there are multiple hydrogens
    neighbors = atoms.neighbors(chiral_atom)
    hydrogens = [atom for atom in neighbors if atom.element == "H"]
    if len(hydrogens) < 3:
        tree = nx.bfs_tree(atoms, chiral_atom)
        # Generate the list of shortest paths in the molecule, neglecting the trivial path [chiral_atom]
        paths = sorted(nx.single_source_shortest_path(tree, chiral_atom, depth).values(), reverse=True, key=lambda x: map(lambda at: at.mass, x))[:-1]
        while paths:
            # Pop the first element (highest priority path) from the list of paths and remove any duplicates.
            path = paths.pop(0)
            # print "Path considered:", path
            paths_no_dups = [unpruned for unpruned in paths if unpruned != path]
            # print "Paths:", paths
            # print "Paths without dups:", paths_no_dups
            # If there are any duplicates, the paths list will be smaller and we can't resolve a highest priority
            if len(paths_no_dups) != len(paths):
                paths = paths_no_dups
            # Otherwise, the path is higher priority than all the other paths, so its second atom is the neighbour with
            # highest priority.
            else:
                # print "Best path:", path
                ranked_atom = path[1]
                # print "Ranked atom:", ranked_atom
                ordered.append(ranked_atom)
                # Drop all the paths containing our ranked atom.
                paths = [unpruned for unpruned in paths if unpruned[1] is not ranked_atom]
    else:
        ordered = [atom for atom in neighbors if atom.element != "H"]
        ordered += [atom for atom in neighbors if atom.element == "H"]
    return ordered

def write_chirality_file(input_filename, output_filename):
  molecule = ra.parse_topology_file(input_filename)
  atoms = molecule.atoms
  chiral_cands = chiral_candidates(atoms)
  multi_bonds(atoms)
  chiral_centres = {}
  for i, chiral_atom in enumerate(chiral_cands):
    ordered = chiral_order(atoms, chiral_atom)
    if len(ordered) == 4:
      chiral_centres[chiral_atom] = ordered
  with open(output_filename, "w") as output_file:
    for atom in sorted(chiral_centres.keys(), cmp=lambda x, y: cmp(x.index, y.index)):
    # Write out the list of chiral atoms and their CIP-ranked neighbours.
      output_string = "{0:>8d}{1:>8d}{2:>8d}{3:>8d}{4:>8d}\n".format(atom.index + 1, 
                                                                     *[other_atom.index + 1 for other_atom in chiral_centres[atom]])
      output_file.write(output_string)

if __name__ == "__main__":
  write_chirality_file(sys.argv[1], ".chirality_list")
