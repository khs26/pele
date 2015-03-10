import networkx as nx
import exceptions as ex
import itertools as it
import re
import playground.group_rotation.amino_acids as amino
import pele.utils.elements as elem


class AmberTopologyFile(object):
    """
    Python representation of the Amber topology file (.prmtop) described here:
    http://ambermd.org/formats.html

    N.B. The convention in topology files is to count from 1, but Python counts
    from 0. I'm going to count from 0 and do the conversion back by hand where
    necessary. This prevents having to force lists to index from 1.
    """

    def __init__(self, filename=None):
        """
        Create the dictionaries and, if we provide a filename, run the read
        method on that filename.
        """
        self.format_dict = {}
        self.data_dict = {}
        self.ordered_flags = []
        self.exclusion_list = {}
        self.input_filename = None
        if filename:
            self.read_from(filename)

    def read_from(self, filename):
        """
        Read from a .prmtop file to the AmberTopologyFile class.
        """
        self.input_filename = filename
        with open(filename, "r") as topology_file:
            input_lines = []
            for line in topology_file:
                line = line.strip("\n")
                if line.startswith('%VERSION'):
                    continue
                if line.startswith('%FLAG'):
                    input_lines.append([line])
                else:
                    input_lines[-1].append(line)
        for input_line in input_lines:
            # Get rid of %FLAG from front and strip whitespace.
            flag = input_line[0][6:].rstrip(" ")
            # Store the flags in order, for writing in future.
            self.ordered_flags.append(flag)
            # Find length and data type.
            type_string = re.sub(r'%FORMAT\((.*)\)[ ]+',
                                 r'\1',
                                 input_line[1])
            # Use capturing brackets around the regex to include the delimiter.
            data_count, data_type, data_format = re.split(r'([a-zA-Z])',
                                                          type_string)
            # We either have a length or a length and precision (e.g. 10F8.4).
            try:
                data_length, data_precision = re.split(r'\.',
                                                       data_format)
                self.format_dict[flag] = (int(data_count),
                                          data_type,
                                          int(data_length),
                                          int(data_precision))
            except ValueError:
                data_length = data_format
                self.format_dict[flag] = (int(data_count),
                                          data_type,
                                          int(data_length))
            data_string = "".join(input_line[2:])
            # Loop through the data string, taking chunks as long as the data length each time.
            split_string = [data_string[i:i + int(data_length)]
                            for i in range(0, len(data_string), int(data_length))]
            # Store the parsed data in the data_dict.
            self.data_dict[flag] = [self.data_cast(data, data_type) for data in split_string]

    def write_to(self, filename):
        """
        Write to a .prmtop file from the AmberTopologyFile class.
        """
        self.update_amber_pointers()
        with open(filename, 'w') as output:
            for flag in self.ordered_flags:
                # Write the flag
                output.write('%FLAG {}'.format(flag).ljust(80))
                output.write('\n')
                # Write the format string
                format_tup = self.format_dict[flag]
                if len(format_tup) == 3:
                    format_str = "".join([str(i) for i in format_tup])
                else:
                    format_str = "".join([str(i) for i in format_tup[:3]] +
                                         ["."] +
                                         [str(format_tup[3])])
                output.write('%FORMAT({})'.format(format_str).ljust(80))
                output.write('\n')
                # Write the data
                # If we have no data to write, just write a newline (for consistency with
                # topology files from LeaP)
                if len(self.data_dict[flag]) == 0:
                    output.write('\n')
                    continue
                # Build the string
                data_count = format_tup[0]
                regrouped = it.izip_longest(fillvalue=None, *[iter(self.data_dict[flag])] * data_count)
                # Define a string for formatting a single unit of data in each of the format cases
                if format_tup[1].lower() == 'a':
                    # e.g. '{:<4}' for a string of length 4
                    format_one = ''.join(['{:<',
                                          str(format_tup[2]),
                                          '}'])
                elif format_tup[1].lower() == 'i':
                    # e.g. '{:>10d}' for an integer of up to width 10
                    format_one = ''.join(['{:>',
                                          str(format_tup[2]),
                                          'd}'])
                elif format_tup[1].lower() == 'e':
                    # e.g. '{:>15.8e}' for an integer of up to width 10
                    format_one = ''.join(['{:>',
                                          str(format_tup[2]),
                                          '.',
                                          str(format_tup[3]),
                                          'E}'])
                else:
                    raise ValueError
                for line in regrouped:
                    data = filter(lambda x: x is not None, line)
                    output.write((format_one * len(data)).format(*data))
                    output.write('\n')

    @staticmethod
    def data_cast(datum, data_type):
        """
        Casts the datum as defined by type.

        a = string
        i = integer
        e = float (in exponent notation)
        """
        if data_type.lower() == "a":
            return str(datum)
        if data_type.lower() == "i":
            return int(datum)
        if data_type.lower() == "e":
            return float(datum)

    def update_amber_pointers(self):
        """
        Updates self.data_dict["POINTERS"] to reflect any changes to the
        topology file. Also see: http://ambermd.org/prmtop.pdf

        A list of pointers is given below. Those marked with * are not
        currently used. We may want to check that these are set to 0.

        0   NATOM    : total number of atoms
            NTYPES   : total number of distinct atom types
            NBONH    : number of bonds containing hydrogen
            MBONA    : number of bonds not containing hydrogen
            NTHETH   : number of angles containing hydrogen
        5   MTHETA   : number of angles not containing hydrogen
            NPHIH    : number of dihedrals containing hydrogen
            MPHIA    : number of dihedrals not containing hydrogen
            *NHPARM  : currently not used
            *NPARM   : used to determine if addles created prmtop
        10  NNB      : number of excluded atoms
            NRES     : number of residues
            NBONA    : MBONA + number of constraint bonds
            NTHETA   : MTHETA + number of constraint angles
            NPHIA    : MPHIA + number of constraint dihedrals
        15  NUMBND   : number of unique bond types
            NUMANG   : number of unique angle types
            NPTRA    : number of unique dihedral types
            *NATYP   : number of atom types in parameter file, see SOLTY below
            *NPHB    : number of distinct 10-12 hydrogen bond pair types
        20  *IFPERT  : set to 1 if perturbation info is to be read in
            *NBPER   : number of bonds to be perturbed
            *NGPER   : number of angles to be perturbed
            *NDPER   : number of dihedrals to be perturbed
            *MBPER   : number of bonds with atoms completely in perturbed group
        25  *MGPER   : number of angles with atoms completely in perturbed group
            *MDPER   : number of dihedrals with atoms completely in perturbed groups
            IFBOX    : set to 1 if standard periodic box, 2 when truncated octahedral
            NMXRS    : number of atoms in the largest residue
            IFCAP    : set to 1 if the CAP option from edit was specified
        30  NUMEXTRA : number of extra points found in topology
            NCOPY    : number of PIMD slices / number of beads
        """
        # Number of atoms
        self.data_dict["POINTERS"][0] = len(self.data_dict["ATOM_NAME"])
        # Number of atom types
        self.data_dict["POINTERS"][1] = len(set(self.data_dict["ATOM_TYPE_INDEX"]))
        # Number of bonds containing hydrogen
        self.data_dict["POINTERS"][2] = len(self.data_dict["BONDS_INC_HYDROGEN"]) / 3
        # Number of bonds not containing hydrogen
        self.data_dict["POINTERS"][3] = len(self.data_dict["BONDS_WITHOUT_HYDROGEN"]) / 3
        # Number of angles containing hydrogen
        self.data_dict["POINTERS"][4] = len(self.data_dict["ANGLES_INC_HYDROGEN"]) / 4
        # Number of angles not containing hydrogen
        self.data_dict["POINTERS"][5] = len(self.data_dict["ANGLES_WITHOUT_HYDROGEN"]) / 4
        # Number of dihedrals containing hydrogen
        self.data_dict["POINTERS"][6] = len(self.data_dict["DIHEDRALS_INC_HYDROGEN"]) / 5
        # Number of dihedrals not containing hydrogen
        self.data_dict["POINTERS"][7] = len(self.data_dict["DIHEDRALS_WITHOUT_HYDROGEN"]) / 5
        # Number of excluded atoms
        self.data_dict["POINTERS"][10] = sum(self.data_dict["NUMBER_EXCLUDED_ATOMS"])
        # Number of residues
        self.data_dict["POINTERS"][11] = len(self.data_dict["RESIDUE_LABEL"])
        # Number of bonds w/o hydrogen + constraint bonds (constraints aren't implemented in
        # topology files any more).
        self.data_dict["POINTERS"][12] = self.data_dict["POINTERS"][3]
        # Number of angles w/o hydrogen + constraint angles
        self.data_dict["POINTERS"][13] = self.data_dict["POINTERS"][5]
        # Number of dihedrals w/o hydrogen + constraint dihedrals
        self.data_dict["POINTERS"][14] = self.data_dict["POINTERS"][7]
        # Number of bond types
        self.data_dict["POINTERS"][15] = len(self.data_dict["BOND_FORCE_CONSTANT"])
        # Number of angle types
        self.data_dict["POINTERS"][16] = len(self.data_dict["ANGLE_FORCE_CONSTANT"])
        # Number of dihedral types
        self.data_dict["POINTERS"][17] = len(self.data_dict["DIHEDRAL_FORCE_CONSTANT"])
        # Number of SOLTY atom types (this isn't used, but the value seems to be set, I
        # include it for consistency with LeaP).
        self.data_dict["POINTERS"][18] = len(self.data_dict["SOLTY"])
        # Largest number of atoms in any residue
        res_starts = self.data_dict["RESIDUE_POINTER"]
        res_ends = self.data_dict["RESIDUE_POINTER"][1:] + [self.data_dict["POINTERS"][0]]
        self.data_dict["POINTERS"][28] = max([y - x for x, y in zip(res_starts, res_ends)])

    def get_exclusion_list(self):
        if not self.exclusion_list:
            start = 0
            # Loop through the NUMBER_EXCLUDED_ATOMS list. Each number refers to the number of
            # elements to take from the EXCLUDED_ATOMS_LIST.
            #
            # e.g. if the first two elements were 3 and 5 and EXCLUDED_ATOMS_LIST started:
            #        2 4 6 3 4 5 6 7
            #      then atom 1 doesn't interact with the first 3 (atoms 2, 4 and 6) and atom
            #      2 doesn't interact with the next 5 (atoms 3-7).
            #
            # The topology file only records exclusions one way, from lower index to higher index
            # so if atoms 7 and 9 don't interact, then atom 7 has 9 in its exclusion list, but not
            # the other way around.
            #
            # If an atom has no excluded atoms (or only those with a lower index than itself), then
            # NUMBER_EXCLUDED_ATOMS has a 1 and EXCLUDED_ATOMS_LIST has a fictitious atom 0 in the
            # relevant place.
            for atom_id, num_excluded in enumerate(self.data_dict["NUMBER_EXCLUDED_ATOMS"]):
                exclusion_list = self.data_dict["EXCLUDED_ATOMS_LIST"][start:start + num_excluded]
                self.exclusion_list[atom_id] = exclusion_list
                start += num_excluded
        return self.exclusion_list

    def add_exclusions(self, atoms_1, atoms_2):
        """
        Add exclusions between every atom index in atoms_1 and every atom in atoms_2. Update
        self.exclusion_list as appropriate.
        """
        if not self.exclusion_list:
            self.get_exclusion_list()
        # Need to make sure that the indexing is from the lower atom to the higher atom.
        to_add = [(min(atom_1, atom_2), max(atom_1, atom_2))
                  for atom_1 in atoms_1 for atom_2 in atoms_2
                  if atom_2 != atom_1]
        # print "Adding:", to_add
        for atom_pair in to_add:
            # Convert to Python numbering (remember we are counting from zero in Python
            # and this dictionary, but 1 in our atom indices).
            atom_id = atom_pair[0] - 1
            # print "Adding to:", atom_id, self.exclusion_list[atom_id]
            if self.exclusion_list[atom_id] == [0]:
                self.exclusion_list[atom_id] = [atom_pair[1]]
            else:
                self.exclusion_list[atom_id].append(atom_pair[1])
        for atom_id in self.exclusion_list:
            self.exclusion_list[atom_id] = sorted(set(self.exclusion_list[atom_id]))

    def remove_exclusions(self, atoms_1, atoms_2):
        """
        Add exclusions between every atom index in atoms_1 and every atom in atoms_2. Update
        self.exclusion_list as appropriate.
        """
        if not self.exclusion_list:
            self.get_exclusion_list()
        # Need to make sure that the indexing is from the lower atom to the higher atom.
        to_remove = [(min(atom_1, atom_2), max(atom_1, atom_2))
                     for atom_1 in atoms_1 for atom_2 in atoms_2
                     if atom_2 != atom_1]
        print "Removing:", to_remove
        for atom_pair in to_remove:
            # Convert to Python numbering (remember we are counting from zero in Python
            # and this dictionary, but 1 in our atom indices).
            atom_id = atom_pair[0] - 1
            # print "Removing from:", atom_id, self.exclusion_list[atom_id]
            self.exclusion_list[atom_id].remove(atom_pair[1])
            if not self.exclusion_list[atom_id]:
                self.exclusion_list[atom_id] = [0]

    def update_exclusion_list(self):
        """
        Updates self.data_dict["NUMBER_EXCLUDED_ATOMS"] and self.data_dict["EXCLUDED_ATOMS_LIST"]
        to reflect the new exclusion list (if it has changed).

        Remember the topology file counts from 1 and Python (and the internal dictionaries) from 0.
        """
        start = 0
        for atom_id in self.exclusion_list:
            # 0..NATOMS-1
            num_excluded = len(self.exclusion_list[atom_id])
            self.data_dict["NUMBER_EXCLUDED_ATOMS"][atom_id] = num_excluded
            self.data_dict["EXCLUDED_ATOMS_LIST"][start:start+num_excluded] = self.exclusion_list[atom_id]
            start += num_excluded


class Atom(object):
    """ Atom defined from the AMBER topology file. """

    def __init__(self, index, name, mass, amber_atom_type, charge, molecule):
        self.index = index
        self.name = name.strip()
        self.mass = mass
        self.amber_atom_type = amber_atom_type
        self.charge = charge
        self.residue = None
        self.element = self.identify_element(self.mass)
        self.molecule = molecule

    @staticmethod
    def identify_element(mass):
        """ Defines the element, based on the mass. """
        return elem.lookup_element_by_mass(mass)

    def set_residue(self, residue):
        """ Defines the Residue to which the Atom belongs. """
        self.residue = residue

    def __repr__(self):
        return " ".join([str(self.index), self.element, self.name])


class Residue(object):
    """ Residue defined from the AMBER topology file. """

    def __init__(self, index, name, molecule):
        self.index = index
        self.name = name.strip()
        self.molecule = molecule
        self.atoms = None

    def add_atoms(self, atoms):
        """ 
        Adds Atoms in atoms to the Residue and sets the residue for those
        Atoms accordingly.
        """
        self.atoms = atoms
        for atom in atoms:
            atom.set_residue(self)

    def bonded(self, other):
        """
        Tests whether this residue is bonded to another residue.
        """
        return any(nx.edge_boundary(self.molecule.atoms,
                                    self.atoms,
                                    other.atoms)) \
               and not self is other

    def __repr__(self):
        return str(self.index) + " " + self.name


class AmberMolecule(object):
    """ Molecule defined from the AMBER topology file. """

    def __init__(self, topology_data=None, coords=None, coords_file=None):
        self.atoms = nx.Graph()
        self.residues = None
        self.coords = coords
        if topology_data:
            self.init_from_topology(topology_data)
        if coords_file:
            self.read_coords(coords_file)

    def init_from_topology(self, topology_data):
        """
        This function takes a topology data set (topology_data) and returns a Molecule object containing
        a graph of the atoms and residues and their connectivity.
        """
        # Create a list of atoms from the topology data and add them to the molecule's graph.
        # The first element of POINTERS is the atom count.
        atom_props = it.izip(topology_data["ATOM_NAME"],
                             topology_data["MASS"],
                             topology_data["AMBER_ATOM_TYPE"],
                             topology_data["CHARGE"])
        atoms = [Atom(i, *prop, molecule=self) for i, prop in enumerate(atom_props)]
        self.atoms.add_nodes_from(atoms)
        # Create a list of residues and add them to the molecule's list of residues.
        # The 11th element of POINTERS is the residue count.
        residues = [Residue(i, resname, self) for i, resname in enumerate(topology_data["RESIDUE_LABEL"])]
        self.residues = residues
        # Go through the BONDS_INC_HYDROGEN and BONDS_WITHOUT_HYDROGEN lists to extract lists of bonded atoms.
        # Index is i / 3 + 1, because AMBER still uses coordinate indices for runtime speed/backwards compatibility.
        first_atoms = [(x / 3) + 1 for x in topology_data["BONDS_INC_HYDROGEN"][0::3] +
                       topology_data["BONDS_WITHOUT_HYDROGEN"][0::3]]
        second_atoms = [(x / 3) + 1 for x in topology_data["BONDS_INC_HYDROGEN"][1::3] +
                        topology_data["BONDS_WITHOUT_HYDROGEN"][1::3]]
        bond_list = zip(first_atoms, second_atoms)
        # Next put the appropriate atoms into the appropriate residues.  AMBER specifies the first atom index of
        # each residue (i.e. no end, hence the exception).
        starts = [x - 1 for x in topology_data["RESIDUE_POINTER"]]
        ends = starts[1:] + [len(atoms)]
        for i, residue in enumerate(residues):
            residue.add_atoms(atoms[starts[i]:ends[i]])
        # Now go through bond list and create bonds between the relevant atoms.
        for bond in bond_list:
            bonded_atoms = (atoms[(bond[0] - 1)], atoms[(bond[1] - 1)])
            self.atoms.add_edge(*bonded_atoms)

    def read_coords(self, filename):
        field_length = 12
        self.coords = []
        with open(filename, "r") as coords_file:
            # Throw away the first line, since it just contains the name of the molecule.
            coords_file.readline()
            # The next line contains the number of atoms.
            number_of_atoms = int(coords_file.readline())
            # Later lines contain coordinates in 12-character wide fields.
            for line in coords_file:
                line = line.rstrip()
                self.coords += map(float, [line[i:i + field_length] for i in range(0, len(line), field_length)])
            # If the number of coordinates is not equal to 3 * number of atoms, raise a RuntimeError.
            if len(self.coords) != number_of_atoms * 3:
                raise ex.RuntimeError("Number of coordinates in coords file and number of atoms are inconsistent.")


def group_rotation_file(molecule, params, filename):
    """
    Defines group rotation for bonds between atoms defined in params.  This function DOESN'T check
    whether the atoms are actually bonded at any point (so make sure you're referencing appropriate
    atom pairs.

    atoms and residues are the atoms and residues lists returned from create_atoms_and_residues()

    params has the format:

    params[(residue, (atom_x, atom_y))] = (rotation_scale, probability)
    e.g. params[("ARG", ("CA", "CB"))] = (0.5, 0.2) rotates about the ARG CA-CB bond at most 90 degs each 
         time, with probability 0.2.
    """
    # Loop through the parameters, creating GROUPROTATION information for any bonds described by an individual
    # parameter.
    with open(filename, "w") as output:
        for param in params.keys():
            for residue in molecule.residues.nodes():
                print residue.name
                if residue.name == param[0]:
                    rotated_atoms = get_rotated_atoms((residue, param[1]))
                    # Create an identifiable group name consisting of [index]_[res_name]_[atom_1]_[atom_2]
                    group_name = "_".join((str(residue.index),
                                           residue.name.strip(),
                                           rotated_atoms[0].name.strip(),
                                           rotated_atoms[1].name.strip()))
                    output.write(" ".join(("GROUP",
                                           group_name,
                                           str(rotated_atoms[0].index + 1),
                                           str(rotated_atoms[1].index + 1),
                                           str(len(rotated_atoms[2])),
                                           str(params[param][0]),
                                           str(params[param][1]),
                                           "\n")))
                    for atom in rotated_atoms[2]:
                        output.write(str(atom.index + 1) + "\n")
                    output.write("\n")


def group_rotation_dict(molecule, params):
    groups = {}
    for param in params.keys():
        for residue in molecule.residues:
            if residue.name == param[0]:
                rotated_atoms = get_rotated_atoms((residue, param[1]))
                # Create an identifiable group name consisting of [index]_[res_name]_[atom_1]_[atom_2]
                group_name = "_".join((str(residue.index),
                                       residue.name.strip(),
                                       rotated_atoms[0].name.strip(),
                                       rotated_atoms[1].name.strip()))
                groups[group_name] = {}
                groups[group_name]["bond_atom_1"] = rotated_atoms[0].index
                groups[group_name]["bond_atom_2"] = rotated_atoms[1].index
                groups[group_name]["group_atoms"] = [atom.index for atom in rotated_atoms[2]]
                groups[group_name]["max_angle_magnitude"] = params[param][0]
                groups[group_name]["selection_probability"] = params[param][1]
    return groups


def get_rotated_atoms(bond):
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # UPDATE ME TO INCLUDE INTER-RES BONDS
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """
    Returns a list of the atoms in the rotating group for the GROUPROTATION script.
    
    bond has the format:

    bond = (residue, (atom_x, atom_y))
    """
    residue = bond[0]
    for residue_atom in residue.atoms:
        if residue_atom.name == bond[1][0]:
            atom_a = residue_atom
        if residue_atom.name == bond[1][1]:
            atom_b = residue_atom
    # Remove the edge, find the smallest subgraph and then add the edge back
    residue.molecule.atoms.remove_edge(atom_a, atom_b)
    rotating_atoms = sorted(nx.connected_components(residue.molecule.atoms), key=lambda x: len(x))[0]
    residue.molecule.atoms.add_edge(atom_a, atom_b)
    if atom_a in rotating_atoms:
        atom_1 = atom_b
        atom_2 = atom_a
        rotating_atoms.remove(atom_a)
    elif atom_b in rotating_atoms:
        atom_1 = atom_a
        atom_2 = atom_b
        rotating_atoms.remove(atom_b)
    return atom_1, atom_2, rotating_atoms


def default_parameters(topology_filename):
    topology_data = AmberTopologyFile(topology_filename)
    parsed = AmberMolecule(topology_data.data_dict)
    return group_rotation_dict(parsed, amino.def_parameters)


if __name__ == "__main__":
    # import sys
    test_top_file = AmberTopologyFile('../amber/coords.prmtop')
    test_top_file.write_to('dummy')
    # for k, v in topology_data.data_dict.items():
    # print k, ":", v
    # for k, v in topology_data.format_dict.items():
    # print k, ":", v
    parsed = AmberMolecule(test_top_file.data_dict)
    parsed.read_coords('../amber/coords.inpcrd')
    print parsed.coords
    for k, v in test_top_file.get_exclusion_list().items():
        print k+1, sorted(set(v))
    test_top_file.add_exclusions([55, 57], [60, 61, 62])
    for k, v in test_top_file.get_exclusion_list().items():
        print k+1, v
    test_top_file.remove_exclusions([60, 62], [59, 63])
    for k, v in test_top_file.get_exclusion_list().items():
        print k + 1, v
    # for res in parsed.residues:
    #     print res, [res2 for res2 in parsed.residues if res.bonded(res2)]
    # group_rot_dict = group_rotation_dict(parsed, amino.def_parameters)
    #params = default_parameters(sys.argv[1])
    # group_rotation_file(parsed,amino.def_parameters,"atomgroups")
    # print group_rot_dict.keys()
    # print parsed
    # print read_amber_coords(sys.argv[2])
    # for item in group_rot_dict:
    #     print item, group_rot_dict[item]
