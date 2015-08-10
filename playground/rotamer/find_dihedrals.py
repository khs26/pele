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


if __name__ == "__main__":
    import os.path

    topology_data = amber.read_topology(os.path.normpath("/home/khs26/flu.prmtop"))
    molecule = amber.create_molecule(topology_data)
    scs = find_sidechains(molecule, [res for res in molecule.residues.nodes()])
    ress, maps = residue_from_sidechain(scs)
    print ress
    print maps