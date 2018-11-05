import re

import numpy as np

from .Biochemistry import *
from ..cy_modules.cy_gridsearch import cy_find_contacts


def weight_function(value):
    """weight function to score contact distances"""
    return 1.0 / (1.0 + np.exp(5.0 * (value - 4.0)))


def loop_trajectory_grid(s1, s2, indices1, indices2, config, suppl,
                         do_self_interaction=False):
    cutoff, hbondcutoff, hbondcutangle = config
    bonds = suppl[0]
    name_array = suppl[1]

    resid_array = []
    segids = []
    if (do_self_interaction):
        resid_array = suppl[2]
        segids = suppl[3]

    # start = time.time()
    currentFrameContacts = []
    natoms1 = len(s1)
    natoms2 = len(s2)
    # TODO:  why so many conversion?
    pos1 = np.array(np.reshape(s1, (1, natoms1 * 3)), dtype=np.float64)
    pos2 = np.array(np.reshape(s2, (1, natoms2 * 3)), dtype=np.float64)
    xyz1 = np.array(pos1, dtype=np.float32)
    xyz2 = np.array(pos2, dtype=np.float32)
    # 2d array with index of atom1 being the index of the first dimension
    # individual lists contain atom2 indices
    nbList1 = cy_find_contacts(xyz1, natoms1, xyz2, natoms2, cutoff)
    idx1 = 0
    for atom1sNeighbors in nbList1:
        for idx2 in atom1sNeighbors:
            convindex1 = indices1[idx1]  # idx1 converted to global atom indexing
            convindex2 = indices2[idx2]  # idx2 converted to global atom indexing
            # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered,
            # hydrogen bonds can still be detected!
            if re.match("H(.*)", name_array[convindex1]) or re.match("H(.*)", name_array[convindex2]):
                continue
            # check if residues are more than 4 apart, and in the same segment
            if do_self_interaction:
                if (resid_array[convindex1] - resid_array[convindex2]) < 5 and segids[convindex1] == segids[convindex2]:
                    continue
            dvec = pos1[0][3*idx1:3*idx1+3] - pos2[0][3*idx2:3*idx2+3]
            distance = np.sqrt(dvec.dot(dvec))

            weight = weight_function(distance)

            # HydrogenBondAlgorithm
            hydrogenBonds = []
            # FF independent hydrogen bonds
            if (name_array[convindex1][0] in HydrogenBondAtoms.atoms and name_array[convindex2][0] in HydrogenBondAtoms.atoms):
                    # print("hbond? %s - %s" % (type_array[convindex1], type_array[convindex2]))
                    # search for hatom, check numbering in bond!!!!!!!!!!
                    b1 = bonds[convindex1]
                    b2 = bonds[convindex2]

                    # b1 = all_sel[convindex1].bonds
                    # b2 = all_sel[convindex2].bonds
                    # search for hydrogen atoms bound to atom 1
                    bondcount1 = 0
                    hydrogenAtomsBoundToAtom1 = []

                    for b in b1.types:
                        hydrogen = next((x for x in b if x.startswith("H")), 0)
                        # print(b)
                        if hydrogen != 0:
                            # print("h bond to atom1")
                            bondindices1 = b1.to_indices()[bondcount1]
                            # print bondindices1
                            # for j in bondindices1:
                            #     print(type_array[j+1])
                            hydrogenidx = next(
                                (j for j in bondindices1 if name_array[j].startswith("H")), -1)
                            if hydrogenidx != -1:
                                hydrogenAtomsBoundToAtom1.append(hydrogenidx)
                        bondcount1 += 1
                    # search for hydrogen atoms bound to atom 2
                    bondcount2 = 0
                    hydrogenAtomsBoundToAtom2 = []
                    # print(b2)
                    for b in b2.types:
                        hydrogen = next((x for x in b if x.startswith("H")), 0)
                        # print(b)
                        if hydrogen != 0:
                            # print("h bond to atom2")
                            bondindices2 = b2.to_indices()[bondcount2]
                            hydrogenidx = next(
                                (k for k in bondindices2 if name_array[k].startswith("H")), -1)
                            if hydrogenidx != -1:
                                # print(type_array[hydrogenidx])
                                hydrogenAtomsBoundToAtom2.append(hydrogenidx)
                        bondcount2 += 1
                    # check hbond criteria for hydrogen atoms bound to first atom
                    for global_hatom in hydrogenAtomsBoundToAtom1:
                        conv_hatom = np.where(indices1[frame] == global_hatom)[0][0]
                        # print(typeHeavy)
                        #
                        # TODO: FF independent version
                        # if (typeHeavy == AtomHBondType.acc or typeHeavy == AtomHBondType.both) and (distarray[conv_hatom, idx2] <= hbondcutoff):
                        # dist = distarray[conv_hatom, idx2]
                        # dist = np.linalg.norm(sel1.positions[conv_hatom] - sel2.positions[idx2])
                        dist = np.linalg.norm(pos1[0][3*conv_hatom:3*conv_hatom+3] - pos2[0][3*idx2:3*idx2+3])
                        if (dist <= hbondcutoff):
                            donorPosition = s1[idx1]
                            hydrogenPosition = np.array(s1[conv_hatom],
                                                        dtype=np.float64)
                            acceptorPosition = np.array(s2[idx2],
                                                        dtype=np.float64)
                            v1 = hydrogenPosition - acceptorPosition
                            v2 = hydrogenPosition - donorPosition
                            v1norm = np.linalg.norm(v1)
                            v2norm = np.linalg.norm(v2)
                            dot = np.dot(v1, v2)
                            angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                            # print(angle)
                            if angle >= hbondcutangle:
                                # print("new hbond")
                                new_hbond = HydrogenBond(convindex1,
                                                         convindex2,
                                                         global_hatom, dist,
                                                         angle,
                                                         hbondcutoff,
                                                         hbondcutangle)
                                hydrogenBonds.append(new_hbond)

                    for global_hatom in hydrogenAtomsBoundToAtom2:
                        conv_hatom = np.where(indices2[frame] == global_hatom)[0][0]
                        dist = np.linalg.norm(pos1[0][3*idx1:3*idx1+3] - pos2[0][3*conv_hatom:3*conv_hatom+3])
                        if (dist <= hbondcutoff):
                            donorPosition = s2[idx2]
                            hydrogenPosition = np.array(s2[conv_hatom], dtype=np.float64)
                            acceptorPosition = np.array(s1[idx1], dtype=np.float64)
                            v1 = hydrogenPosition - acceptorPosition
                            v2 = hydrogenPosition - donorPosition
                            v1norm = np.linalg.norm(v1)
                            v2norm = np.linalg.norm(v2)
                            dot = np.dot(v1, v2)
                            angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                            if angle >= hbondcutangle:
                                new_hbond = HydrogenBond(convindex2, convindex1, global_hatom, dist, angle,
                                                         hbondcutoff,
                                                         hbondcutangle)
                                hydrogenBonds.append(new_hbond)

            newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1),
                                         int(convindex2),
                                         hydrogenBonds)
            currentFrameContacts.append(newAtomContact)

    return currentFrameContacts
