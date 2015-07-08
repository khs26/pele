"""
-----------
Setup stage
-----------

- Identify residues in the molecule.
- Find dihedrals and store atom identities.
- Choose backbone-dependent, local environment etc. to determine what sort of library data to generate.
- Read in relevant part of rotamer library data: store dictionary for each residue's rotamer states.

----------
Move stage
----------

- Read in coordinates from GMIN.
- Select a set of residues to move.
- Measure rotamer state of residues and their neighbours (if appropriate).
- Select a new conformation from the rotamer library, or a random configuration, with probability given by
  Good-Turing frequency estimation.
- Write coordinates for GMIN.

------
Future
------

- Decomposed energies for dihedrals wrt other sidechain dihedrals.

"""

