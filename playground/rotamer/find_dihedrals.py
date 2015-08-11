import pele.amber.read_amber as amber
import networkx as nx
import playground.group_rotation.chirality as chirality
from playground.rotamer.measure_dihedral import dihedral_with_symmetry

class Dihedral(object):
    def __init__(self, dihedral_atoms):
        self.atoms = dihedral_atoms
        self.residue = dihedral_atoms[0].residue
        self.atom_graph = dihedral_atoms[0].molecule.atoms
        self.get_moving_atoms()
    def measure_dihedral(self, coords):
        return dihedral_with_symmetry(coords, self)
    def change_dihedral(self, start_coords):
        current_angle = self.measure_dihedral(start_coords)
    def get_moving_atoms(self):
        """ Returns a list of moving atoms for a given dihedral.

        :param dihedral: List of atoms which define the dihedral.
        :param atom_graph: Atom graph.
        :return moving_atoms: list of moving atoms.
        """
        self.atom_graph.remove_edge(dihedral[1], dihedral[2])
        moving_atoms = [g for g in nx.connected_components(self.atom_graph) if dihedral[2] in g][0]
        self.moving_atoms = [atom for atom in moving_atoms if not isinstance(atom, chirality.GhostAtom)]
        self.atom_graph.add_edge(dihedral[1], dihedral[2])
        return self.moving_atoms


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
    from playground.rotamer.measure_dihedral import dihedrals_with_symmetry
    import cProfile

    pr = cProfile.Profile()
    pr.enable()
    topology_data = amber.read_topology(os.path.normpath("/home/khs26/flu.prmtop"))
    coords = np.array(amber.read_amber_coords(os.path.normpath("/home/khs26/flu.inpcrd"))).reshape((-1, 3))
    molecule = amber.create_molecule(topology_data)
    ress, maps = molecule.identify_residues()
    dihedrals = []
    for v in map_dihedrals(ress, maps).values():
        for dihedral in v:
            dihedrals.append(Dihedral(dihedral))
    for dihe in dihedrals:
        print dihe.residue, dihe.atoms, dihe.moving_atoms
        print "Angle:", dihe.measure_dihedral(coords)
    pr.disable()
    pr.print_stats('cumulative')