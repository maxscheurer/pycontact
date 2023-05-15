import numpy as np
from PyContact.cy_modules.cy_gridsearch import cy_find_contacts
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist


def test_find_contacts():
    # pass positions and distance to cython
    natoms1 = 1800
    natoms2 = 2100
    box_dimensions = [50, 55, 45]
    # generate random positions in a box of size box_dimensions
    np.random.seed(42)
    pos1 = np.random.rand(natoms1, 3) * box_dimensions
    pos2 = np.random.rand(natoms2, 3) * box_dimensions
    xyz1 = np.array(np.reshape(pos1, (1, natoms1 * 3)), dtype=np.float32)
    xyz2 = np.array(np.reshape(pos2, (1, natoms2 * 3)), dtype=np.float32)
    
    cutoff = 5.0
    # 2d array with index of atom1 being the index of the first dimension
    # individual lists contain atom2 indices
    nbList1 = cy_find_contacts(xyz1, natoms1, xyz2, natoms2, cutoff)
    
    tree = KDTree(pos2)
    # Query the KD-Tree to find pairs of points within the cutoff distance
    pairs = tree.query_ball_point(pos1, r=cutoff)

    # compute pair list using scipy.spatial distance matrix
    dist_mat = cdist(pos1, pos2)
    # get indices of atoms that are within cutoff
    pair_list = np.where(dist_mat <= cutoff)
    assert len(pair_list[0]) > 0
    assert len(pair_list[1]) > 0
    # check if the number of contacts is the same
    for (p1, p2) in zip(pair_list[0], pair_list[1]):
        assert p2 in nbList1[p1]
        assert p2 in pairs[p1]
    
    for idx, nb_list in enumerate(nbList1[:natoms1]):
        assert sorted(nb_list) == sorted(pairs[idx])