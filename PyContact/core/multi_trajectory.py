from __future__ import print_function
import os
import re
import time
from copy import deepcopy
import itertools

import MDAnalysis
from MDAnalysis.analysis import distances
import numpy as np

from .Biochemistry import *
from .LogPool import *
from ..cy_modules.cy_gridsearch import cy_find_contacts

def weight_function(value):
    """weight function to score contact distances"""
    return 1.0 / (1.0 + np.exp(5.0 * (value - 4.0)))


def chunks(seq, num):
    """splits the list seq in num (almost) equally sized chunks."""
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


class ConvBond(object):
    """Python object of MDAnalysis bond for running jobs in parallel."""
    def __init__(self, bonds):
        super(ConvBond, self).__init__()
        self.types = []
        self.indices = []
        for b in bonds:
            self.indices.append(deepcopy(b.indices))
            self.types.append(deepcopy(b.type))

    def types(self):
        return self.types

    def to_indices(self):
        return self.indices

# newer version with gridsearch
def loop_trajectory_grid(sel1c, sel2c, indices1, indices2, config, suppl, selfInteraction):
    cutoff, hbondcutoff, hbondcutangle = config
    # resname_array = comm.bcast(resname_array, root=0)
    # resid_array = comm.bcast(resid_array, root=0)
    # name_array = comm.bcast(name_array, root=0)
    bonds = suppl[0]
    # segids = comm.bcast(segids, root=0)
    # backbone = comm.bcast(backbone, root=0)
    name_array = suppl[1]

    resid_array = []
    segids = []
    if (selfInteraction):
        resid_array = suppl[2]
        segids = suppl[3]

    allRankContacts = []
    # start = time.time()
    for s1, s2 in zip(sel1c, sel2c):
        frame = 0
        currentFrameContacts = []
        natoms1 = len(s1)
        natoms2 = len(s2)
        pos1 = np.reshape(s1, (1, natoms1 * 3))
        pos2 = np.reshape(s2, (1, natoms2 * 3))
        xyz1 = np.array(pos1, dtype=np.float32)
        xyz2 = np.array(pos2, dtype=np.float32)
        # 2d array with index of atom1 being the index of the first dimension
        # individual lists contain atom2 indices
        nbList1 = cy_find_contacts(xyz1, natoms1, xyz2, natoms2, cutoff)
        # nbList1 = res[:natoms1]
        # nbList2 = res[natoms1:]
        idx1 = 0
        for atom1sNeighbors in nbList1:
            for idx2 in atom1sNeighbors:
                convindex1 = indices1[frame][idx1]  # idx1 converted to global atom indexing
                convindex2 = indices2[frame][idx2]  # idx2 converted to global atom indexing
                # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered,
                # hydrogen bonds can still be detected!
                if re.match("H(.*)", name_array[convindex1]) or re.match("H(.*)", name_array[convindex2]):
                    continue
                    # distance between atom1 and atom2
                # check if residues are more than 4 apart, and in the same segment
                if selfInteraction:
                    if (resid_array[convindex1] - resid_array[convindex2]) < 5 and segids[convindex1] == segids[convindex2]:
                        continue
                # distance = distarray[idx1, idx2]
                # weight = weight_function(distance)
                distance = np.linalg.norm(pos1[0][3*idx1:3*idx1+3] - pos2[0][3*idx2:3*idx2+3])
                # if (distance - distarray[idx1, idx2]) > 0.001:
                #     print("Error in distance calculations!")
                #     return
                # # print(convindex1, convindex2, distance, distarray[idx1, idx2])
                # if (distance > cutoff):
                #     print("Distances must be smaller/equal cutoff!")
                #     return
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

                        # new code
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
                                    # print(type_array[hydrogenidx])
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
                                hydrogenPosition = s1[conv_hatom]
                                acceptorPosition = s2[idx2]
                                v1 = hydrogenPosition - acceptorPosition
                                v2 = hydrogenPosition - donorPosition
                                v1norm = np.linalg.norm(v1)
                                v2norm = np.linalg.norm(v2)
                                dot = np.dot(v1, v2)
                                angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                                # print(angle)
                                if angle >= hbondcutangle:
                                    # print("new hbond")
                                    new_hbond = HydrogenBond(convindex1, convindex2, global_hatom, dist, angle,
                                                             hbondcutoff,
                                                             hbondcutangle)
                                    hydrogenBonds.append(new_hbond)
                                    # print(str(convindex1) + " " + str(convindex2)
                                    # print("hbond found: %d,%d,%d"%(convindex1,global_hatom,convindex2))
                                    # print(angle)
                        for global_hatom in hydrogenAtomsBoundToAtom2:
                            conv_hatom = np.where(indices2[frame] == global_hatom)[0][0]
                            # TODO: FF independent version
                            # if (typeHeavy == AtomHBondType.acc or typeHeavy == AtomHBondType.both) and (distarray[idx1, conv_hatom] <= hbondcutoff):
                            # FIXME: WTF?
                            # if (distarray[conv_hatom, idx2] <= hbondcutoff):
                            # dist = distarray[idx1, conv_hatom]
                            # dist = np.linalg.norm(sel1.positions[idx1] - sel2.positions[conv_hatom])
                            dist = np.linalg.norm(pos1[0][3*idx1:3*idx1+3] - pos2[0][3*conv_hatom:3*conv_hatom+3])
                            if (dist <= hbondcutoff):
                                donorPosition = s2[idx2]
                                hydrogenPosition = s2[conv_hatom]
                                acceptorPosition = s1[idx1]
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
                                    # print str(convindex1) + " " + str(convindex2)
                                    # print "hbond found: %d,%d,%d"%(convindex2,global_hatom,convindex1)
                                    # print angle
                                    # finalize
                newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1),
                                             int(convindex2),
                                             hydrogenBonds)
                currentFrameContacts.append(newAtomContact)
            idx1 += 1
        allRankContacts.append(currentFrameContacts)
        frame += 1

    return allRankContacts

def loop_trajectory(sel1c, sel2c, indices1, indices2, config, suppl, selfInteraction):
    """Invoked to analyze trajectory chunk for contacts as a single thread."""
    # print(len(sel1c) , len(sel2c))
    # indices1 = suppl[0]
    # indices2 = suppl[1]
    cutoff, hbondcutoff, hbondcutangle = config
    # resname_array = comm.bcast(resname_array, root=0)
    # resid_array = comm.bcast(resid_array, root=0)
    # name_array = comm.bcast(name_array, root=0)
    bonds = suppl[0]
    # segids = comm.bcast(segids, root=0)
    # backbone = comm.bcast(backbone, root=0)
    name_array = suppl[1]

    resid_array = []
    segids = []
    if (selfInteraction):
        resid_array = suppl[2]
        segids = suppl[3]

    allRankContacts = []
    # start = time.time()
    for s1, s2 in zip(sel1c, sel2c):
        frame = 0
        currentFrameContacts = []
        result = np.ndarray(shape=(len(s1), len(s2)), dtype=float)
        distarray = distances.distance_array(s1, s2, box=None, result=result)
        contacts = np.where(distarray <= cutoff)
        for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
            convindex1 = indices1[frame][idx1]  # idx1 converted to global atom indexing
            convindex2 = indices2[frame][idx2]  # idx2 converted to global atom indexing
            # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered,
            # hydrogen bonds can still be detected!
            if re.match("H(.*)", name_array[convindex1]) or re.match("H(.*)", name_array[convindex2]):
                continue
            if selfInteraction:
                if (resid_array[convindex1] - resid_array[convindex2]) < 5 and segids[convindex1] == segids[convindex2]:
                    continue
            # distance between atom1 and atom2
            distance = distarray[idx1, idx2]
            weight = weight_function(distance)
            hydrogenBonds = []
            if (name_array[convindex1][0] in HydrogenBondAtoms.atoms and name_array[convindex2][0] in HydrogenBondAtoms.atoms):
                    # print("hbond? %s - %s" % (type_array[convindex1], type_array[convindex2]))
                    # search for hatom, check numbering in bond!!!!!!!!!!
                    b1 = bonds[convindex1]
                    b2 = bonds[convindex2]
                    bondcount1 = 0
                    hydrogenAtomsBoundToAtom1 = []
                    # new code
                    for b in b1.types:
                        # b = bnd.type
                        hydrogen = next((xx for xx in b if xx.startswith("H")), 0)
                        # print(b)
                        if hydrogen != 0:
                            # print("h bond to atom1")
                            bondindices1 = b1.to_indices()[bondcount1]
                            # print bondindices1
                            # for j in bondindices1:
                            #     print(self.type_array[j+1])
                            hydrogenidx = next(
                                (j for j in bondindices1 if name_array[j].startswith("H")), -1)
                            if hydrogenidx != -1:
                                # print(self.type_array[hydrogenidx])
                                hydrogenAtomsBoundToAtom1.append(hydrogenidx)
                        bondcount1 += 1
                    # search for hydrogen atoms bound to atom 2
                    bondcount2 = 0
                    hydrogenAtomsBoundToAtom2 = []
                    for b in b2.types:
                        # b = bnd2.type
                        hydrogen = next((xx for xx in b if xx.startswith("H")), 0)
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

                    for global_hatom in hydrogenAtomsBoundToAtom1:
                        conv_hatom = np.where(indices1[frame] == global_hatom)[0][0]
                        dist = distarray[conv_hatom, idx2]
                        if dist <= hbondcutoff:
                            donorPosition = s1[idx1]
                            hydrogenPosition = s1[conv_hatom]
                            acceptorPosition = s2[idx2]
                            v1 = hydrogenPosition - acceptorPosition
                            v2 = hydrogenPosition - donorPosition
                            v1norm = np.linalg.norm(v1)
                            v2norm = np.linalg.norm(v2)
                            dot = np.dot(v1, v2)
                            angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                            if angle >= hbondcutangle:
                                new_hbond = HydrogenBond(convindex1, convindex2, global_hatom, dist, angle,
                                                         hbondcutoff,
                                                         hbondcutangle)
                                hydrogenBonds.append(new_hbond)
                            # print str(convindex1) + " " + str(convindex2)
                            # print "hbond found: %d,%d,%d"%(convindex1,global_hatom,convindex2)
                            # print angle
                    for global_hatom in hydrogenAtomsBoundToAtom2:
                        conv_hatom = np.where(indices2[frame] == global_hatom)[0][0]
                        dist = distarray[idx1, conv_hatom]
                        if dist <= hbondcutoff:
                            donorPosition = s2[idx2]
                            hydrogenPosition = s2[conv_hatom]
                            acceptorPosition = s1[idx1]
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
            newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1), int(convindex2),
                                         hydrogenBonds)
            currentFrameContacts.append(newAtomContact)
        allRankContacts.append(currentFrameContacts)
        frame += 1
    return allRankContacts


def run_load_parallel(nproc, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
    """Invokes nproc threads to run trajectory loading and contact analysis in parallel."""
    # nproc = int(self.settingsView.coreBox.value())
    pool = LoggingPool(nproc)
    # manager = multiprocessing.Manager()
    # d=manager.list(trajArgs)

    # load psf and dcd
    u = MDAnalysis.Universe(psf, dcd)
    # define selections according to sel1text and sel2text

    selfInteraction = False
    if sel2text == "self":
        sel1 = u.select_atoms(sel1text)
        sel2 = u.select_atoms(sel1text)
        selfInteraction = True
    else:
        sel1 = u.select_atoms(sel1text)
        sel2 = u.select_atoms(sel2text)

    # write properties of all atoms to lists
    all_sel = u.select_atoms("all")
    backbone_sel = u.select_atoms("backbone")
    resname_array = []
    resid_array = []
    name_array = []
    bonds = []
    segids = []
    backbone = []
    for atom in all_sel.atoms:
        resname_array.append(atom.resname)
        resid_array.append(atom.resid)
        name_array.append(atom.name)
        bonds.append(ConvBond(atom.bonds))
        segids.append(atom.segid)
    for atom in backbone_sel:
        backbone.append(atom.index)

    if (len(sel1.atoms) == 0 or len(sel2.atoms) == 0):
        raise Exception

    sel1coords = []
    sel2coords = []
    start = time.time()
    indices1 = []
    indices2 = []

    for ts in u.trajectory:
        # define selections according to sel1text and sel2text
        if "around" in sel1text:
            sel1 = u.select_atoms(sel1text)
        if "around" in sel2text:
            sel2 = u.select_atoms(sel2text)
        # write atomindices for each selection to list
        sel1coords.append(sel1.positions)
        sel2coords.append(sel2.positions)
        # tempindices1 = []
        # for at in sel1.atoms:
        #     tempindices1.append(at.index)
        # tempindices2 = []
        # for at in sel2.atoms:
        #     tempindices2.append(at.index)
        indices1.append(sel1.indices)
        indices2.append(sel2.indices)

    # contactResults = []
    # loop over trajectory
    # totalFrameNumber = len(u.trajectory)
    # start = time.time()
    sel1c = chunks(sel1coords, nproc)
    sel2c = chunks(sel2coords, nproc)
    sel1ind = chunks(indices1, nproc)
    sel2ind = chunks(indices2, nproc)
    # print(len(sel1ind), len(sel2ind))
    # show trajectory information and selection information
    # print("trajectory with %d frames loaded" % len(u.trajectory))

    print("Running on %d cores" % nproc)
    results = []
    rank = 0
    for c in zip(sel1c, sel2c, sel1ind, sel2ind):
        if (selfInteraction):
            results.append(pool.apply_async(loop_trajectory_grid, args=(c[0], c[1], c[2], c[3],
                                                               [cutoff, hbondcutoff, hbondcutangle],
                                                               [bonds, name_array, resid_array, segids], selfInteraction)))
        else:
            results.append(pool.apply_async(loop_trajectory_grid, args=(c[0], c[1], c[2], c[3],
                                                               [cutoff, hbondcutoff, hbondcutangle],
                                                               [bonds, name_array], selfInteraction)))
        rank += 1
    pool.close()
    pool.join()
    # stop = time.time()
    allContacts = []
    for res in results:
        rn = res.get()
        # print(len(rn))
        allContacts.extend(rn)
    # pickle.dump(allContacts,open("parallel_results.dat","w"))
    # print("frames: ", len(allContacts))
    # print("time: ", str(stop-start), rank)
    return [allContacts, resname_array, resid_array, name_array, segids, backbone]
