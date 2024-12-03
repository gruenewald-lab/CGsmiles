import numpy as np
from numpy import unravel_index
from scipy.spatial.distance import pdist, squareform

def u_vect(vect):
    """
    Compute unit vector of vector.

    Parameters
    ----------
    vect: np.array

    Returns
    -----------
    np.array
    """
    u_vect = vect/np.linalg.norm(vect)
    return u_vect

def vector_angle_degrees(v1, v2):
    """
    Compute the angle between two vectors
    in degrees and between 0 180 degrees.

    Parameters
    ----------
    v1: np.array
    v2: np.array

    Returns
    ---------
    float
      angle in degrees
    """
    angle = np.degrees(np.arccos(np.clip(np.dot(u_vect(v1), u_vect(v2)),-1, 1)))
    return angle

def dihedral_angle_between(a1, a2, a3, a4):
    """
    Compute dihedral angle between positions.
    """
    r1 = a1 - a2
    r2 = a3 - a2
    r3 = a3 - a4
    cross1 = np.cross(r1, r2)
    cross2 = np.cross(r2, r3)
    n1 = u_vect(cross1)
    n2 = u_vect(cross2)
    dih = vector_angle_degrees(n1, n2)
    return dih

def rotate(positions, angle, origin=np.array([0, 0])):
    """
    Rotate positions by angle in 2D about the origin.
    The angle is given in radians.
    """
    positions = positions - origin
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])
    rotated_positions = np.dot(positions, rotation_matrix.T)
    rotated_positions = rotated_positions + origin
    return rotated_positions

def rotate_degrees(position, angle, origin=np.array([0,0])):
    """
    Wrapper for rotating about the origin but with angle
    given in degrees.
    """
    angle = np.deg2rad(angle)
    return rotate(position, angle, origin)

def rotate_to_axis(positions, align_with, diagonal=None):
    """
    Rotate position array to aling with one of the principle axis
    or the diagonal.
    """
    # get all the distances
    distances = squareform(pdist(positions, 'euclidean'))
    # find pair belonging to largest distance
    pair = unravel_index(distances.argmax(), distances.shape)
    # get vector
    vector = positions[pair[0]] - positions[pair[1]]
    # compute angle
    if align_with == 'diag':
        diagonal = np.array(diagonal)
    elif align_with == 'x':
        diagonal = np.array([1,0])
    elif align_with == 'y':
        diagonal = np.array([0, 1])

    angle = np.arctan2(diagonal[1], diagonal[0]) - np.arctan2(vector[1], vector[0])
    rotated_positions = rotate(positions, angle)
    return rotated_positions
