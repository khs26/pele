import numpy as np

def dihedral_angle(coords, dihedral):
    """
    Measure dihedral angle for this set of atoms.
    :param coords:
    :param dihedral:
    :return angle: measured angle
    """
    atom_coords = [coords[x.index] for x in dihedral]
    b1, b2, b3 = [atom_coords[i+1] - atom_coords[i] for i in range(0, 3)]
    b2_b3 = np.cross(b2, b3)
    b1_b2 = np.cross(b1, b2)
    angle = np.arctan2(np.linalg.norm(b2) * np.dot(b1, b2_b3), np.dot(b1_b2, b2_b3))
    return angle