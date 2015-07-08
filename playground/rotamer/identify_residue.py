import pele.amber.read_amber as amber
import playground.group_rotation.chirality as chir
import numpy as np
import networkx as nx


def res_atom_graph(molecule_graph, residues):
    """returns the graph of atoms in residues

    :param molecule: molecule containing residue and atom graphs
    :type molecule: amber.Molecule
    :param residues: residues
    :type residues: amber.Residue
    :return: graph containing only atoms in residues
    :rtype: nx.Graph
    """
    res_atoms = []
    for res in residues:
        res_atoms += res.atoms
    return molecule_graph.atoms.subgraph(res_atoms)


def find_sidechains(molecule_graph, residues):
    # Identify chiral atoms
    atoms = molecule.atoms
    chiral_centres = chir.get_chiral_sets(atoms)
    # Identify sidechains (Ca-Cb-X), apart from proline and glycine.
    sidechains = {}
    for k, v in chiral_centres.items():
        carbons = [atom for atom in v if atom.element == 'C']
        amides = [carbon for carbon in carbons
                  if any([type(nb) == chir.GhostAtom and nb.element == 'O' for nb in nx.neighbors(atoms, carbon)])
                  and any([nb.element == 'N' or nb.element == 'O' for nb in nx.neighbors(atoms, carbon)])]
        nbs_n = [nb for nb in v if nb.element == 'N']
        if amides and nbs_n:
            amide_bond = (k, amides[0])
            n_bond = (k, nbs_n[0])
            h_bond = (k, [h for h in nx.neighbors(atoms, k) if h.element == 'H'][0])
            # Now find sidechains by cutting the Ca-C, Ca-N and Ca-H bonds
            atoms.remove_edges_from([amide_bond, n_bond, h_bond])
            sidechain_atoms = [atom for atom in [comp for comp in nx.connected_components(atoms) if k in comp][0]
                               if type(atom) != chir.GhostAtom]
            atoms.add_edges_from([amide_bond, n_bond, h_bond])
            if not any([k in cycle for cycle in nx.cycle_basis(atoms.subgraph(sidechain_atoms))]):
                sidechains[k] = atoms.subgraph(sidechain_atoms)
    return sidechains

def residue_from_sidechain():
    pass

if __name__ == "__main__":
    topology_data = amber.read_topology("/home/khs26/coords.prmtop")
    molecule = amber.create_molecule(topology_data)
    scs = find_sidechains(molecule, [res for res in molecule.residues.nodes() if 'T' in res.name])
    for sc in scs:
        for sc2 in scs:
            if sc == sc2:
                continue
            if nx.is_isomorphic(scs[sc], scs[sc2], node_match=(lambda x, y: x['element'] == y['element'])):
                print sc, sc2, "are isomorphic"
                print sc.residue, sc2.residue
                print nx.neighbors(scs[sc], sc)
    # test_coords = np.array(amber.read_amber_coords("/home/khs26/coords.inpcrd"))