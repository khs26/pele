import pele.amber.read_amber as amber
import playground.group_rotation.chirality as chir
import residue_sidechains as res_scs
from identify_residue import find_sidechains, residue_from_sidechain
import networkx as nx


def expand_dihedral(atom1, atom2):
    """ Finds the appropriate neighbouring atoms for atom1 and atom2, to define the dihedral properly.

    This first tries to find the atom with the highest rank according to the Cahn-Ingold-Prelog rules for determining
    chirality. If there is a tie (e.g. valine, leucine, tyrosine), it instead takes one of the equivalent atoms and
    notes the symmetry of the atoms.

    :param atom1: atom 1 of pair
    :param atom2: atom 2 of pair
    :return:
    """
    atoms = atom1.molecule.atoms
    sym0, sym3 = False, False
    # Try the chirality test
    try:
        ordered = [i for i in chir.chiral_order(atoms, atom1) if i != atom2]
        atom0 = ordered[0]
        if len(ordered) < 3:
            sym0 = True
    except IndexError:
        atom0 = sorted(nx.neighbors(atoms, atom1), key=lambda x: x.mass)[-1]
    try:
        ordered = [i for i in chir.chiral_order(atoms, atom2) if i != atom1]
        atom3 = ordered[0]
        if "VAL" in atom3.residue.name:
            print "Ordered:", ordered
            exit()
        if len(ordered) < 3:
            sym3 = True
    except IndexError:
        atom3 = sorted(nx.neighbors(atoms, atom2), key=lambda x: x.mass)[-1]
    return atom0, atom1, atom2, atom3, sym0, sym3


def map_dihedrals(residue_identities, residue_atom_map):
    """ Generates tuples of atoms in each dihedral for every residue.

    :param residue_identities: Dictionary of atom identities.
    :param residue_atom_map: Dictionary of atom mappings from the molecule.atoms graph to corresponding
    residue_sidechains atom names.
    :return dihedral_dict: Dictionary of dihedrals for each residue in molecule.residues. Entries are lists of lists,
    where each element of the sublist is an atom in molecule.residue.atoms.
    """
    from residue_sidechains import dihedral_chain
    dihedral_dict = {}
    for residue, res_id in residue_identities.items():
        # If the residue doesn't have an entry in dihedral_chain, skip to the next residue.
        if not dihedral_chain[res_id]:
            print "Skipping:", residue, res_id
            continue
        unmapped_dihedrals = [dihedral_chain[res_id][i:i+4] for i in range(0, len(dihedral_chain[res_id]) - 3)]
        mapped_dihedrals = []
        for dihedral in unmapped_dihedrals:
            # We need to convert "N-1" to the N neighbouring "C0". The others can be converted using the map.
            atom_map = residue_atom_map[residue]
            mapped_dihedral = [[nb for nb in nx.neighbors(residue.molecule.atoms, atom_map["C0"]) if nb.element == "N"][0]
                               if atom_name == "N-1" else atom_map[atom_name] for atom_name in dihedral]
            mapped_dihedrals.append(mapped_dihedral)
        dihedral_dict[residue] = mapped_dihedrals
    return dihedral_dict


if __name__ == "__main__":
    import os.path
    import numpy as np

    topology_data = amber.read_topology(os.path.normpath("/home/khs26/flu.prmtop"))
    coords = np.array(amber.read_amber_coords(os.path.normpath("/home/khs26/flu.inpcrd"))).reshape((-1, 3))
    molecule = amber.create_molecule(topology_data)
    scs = find_sidechains(molecule)
    print scs
    ress, maps = residue_from_sidechain(scs)
    print ress
    print maps
    for k, v in map_dihedrals(ress, maps).items():
        print k, [map(lambda x: (x.name, x.index, coords[x.index]), dihedral) for dihedral in v]