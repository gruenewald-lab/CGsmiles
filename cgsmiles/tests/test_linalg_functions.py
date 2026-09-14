import numpy as np
import pytest
from cgsmiles.linalg_functions import (u_vect,
                                       vector_angle_degrees,
                                       dihedral_angle_between,
                                       rotate,
                                       rotate_degrees,
                                       rotate_to_axis)


@pytest.mark.parametrize('vect', (
                         np.array([3.0, 4.0]),
                         np.array([1.0, 0.0, 0.0]),
                         np.array([-2.0, 5.0, -7.0]),
))
def test_u_vect_has_unit_norm(vect):
    unit = u_vect(vect)
    assert np.isclose(np.linalg.norm(unit), 1.0)
    assert np.allclose(unit, vect / np.linalg.norm(vect))


def test_u_vect_preserves_direction():
    vect = np.array([3.0, 4.0])
    unit = u_vect(vect)
    assert np.allclose(unit, [0.6, 0.8])


@pytest.mark.parametrize('v1, v2, expected', (
                         (np.array([1.0, 0.0]), np.array([1.0, 0.0]), 0.0),
                         (np.array([1.0, 0.0]), np.array([-1.0, 0.0]), 180.0),
                         (np.array([1.0, 0.0]), np.array([0.0, 1.0]), 90.0),
                         (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), 90.0),
                         # scaling a vector must not change the angle
                         (np.array([2.0, 0.0]), np.array([0.0, 5.0]), 90.0),
))
def test_vector_angle_degrees(v1, v2, expected):
    assert np.isclose(vector_angle_degrees(v1, v2), expected)


def test_vector_angle_degrees_is_clipped_for_parallel_vectors():
    # identical vectors give a dot product of exactly 1 for the unit
    # vectors; this must not raise from arccos domain errors
    v1 = np.array([1.0, 2.0, 3.0])
    assert np.isclose(vector_angle_degrees(v1, v1), 0.0)
    assert np.isclose(vector_angle_degrees(v1, -v1), 180.0)


@pytest.mark.parametrize('a4, expected', (
                         (np.array([-1.0, 0.0, 1.5]), 180.0),
                         (np.array([1.0, 0.0, 1.5]), 0.0),
                         (np.array([0.0, 1.0, 1.5]), 90.0),
))
def test_dihedral_angle_between(a4, expected):
    a1 = np.array([1.0, 0.0, -0.5])
    a2 = np.array([0.0, 0.0, 0.0])
    a3 = np.array([0.0, 0.0, 1.0])
    assert np.isclose(dihedral_angle_between(a1, a2, a3, a4), expected)


def test_rotate_quarter_turn_about_origin():
    positions = np.array([1.0, 0.0])
    rotated = rotate(positions, np.pi / 2)
    assert np.allclose(rotated, [0.0, 1.0], atol=1e-8)


def test_rotate_full_turn_is_identity():
    positions = np.array([[1.0, 0.0], [0.0, 2.0], [-3.0, 1.0]])
    rotated = rotate(positions, 2 * np.pi)
    assert np.allclose(rotated, positions, atol=1e-8)


def test_rotate_batch_of_positions():
    positions = np.array([[1.0, 0.0], [0.0, 1.0]])
    rotated = rotate(positions, np.pi / 2)
    assert np.allclose(rotated, [[0.0, 1.0], [-1.0, 0.0]], atol=1e-8)


def test_rotate_about_non_origin_point():
    positions = np.array([1.0, 0.0])
    rotated = rotate(positions, np.pi / 2, origin=np.array([1.0, 0.0]))
    assert np.allclose(rotated, [1.0, 0.0], atol=1e-8)


def test_rotate_degrees_matches_rotate_in_radians():
    positions = np.array([[2.0, 1.0], [-1.0, 3.0]])
    by_degrees = rotate_degrees(positions, 90)
    by_radians = rotate(positions, np.pi / 2)
    assert np.allclose(by_degrees, by_radians)


def test_rotate_to_axis_aligns_longest_distance_pair():
    positions = np.array([[0.0, 0.0], [1.0, 0.5], [2.0, 0.0], [-1.0, -0.5]])
    rotated = rotate_to_axis(positions, np.array([0.0, 1.0]))

    # the two points that were furthest apart before rotation must now
    # share the same x-coordinate, i.e. lie along the y-axis
    distances = np.linalg.norm(positions[:, None] - positions[None, :], axis=-1)
    i, j = np.unravel_index(distances.argmax(), distances.shape)
    assert np.isclose(rotated[i][0], rotated[j][0])
    # and the overall pairwise distance is preserved by the rotation
    assert np.isclose(np.linalg.norm(rotated[i] - rotated[j]),
                      np.linalg.norm(positions[i] - positions[j]))


def test_rotate_to_axis_with_already_aligned_positions_stays_on_axis():
    # points already spread out along the x-axis, aligned with x-axis target;
    # the rotation may still flip the sign (180 degree ambiguity) since
    # alignment only cares about the line, not the direction along it
    positions = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    rotated = rotate_to_axis(positions, np.array([1.0, 0.0]))
    assert np.allclose(rotated[:, 1], 0.0, atol=1e-8)
