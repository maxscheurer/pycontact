from __future__ import print_function
import os
import re
import sys, time
from copy import deepcopy
import itertools

import MDAnalysis
from MDAnalysis.analysis import distances

from .Biochemistry import *
from .LogPool import *

def weight_function(value):
    return (1.0) / (1.0 + np.exp(5.0 * (value - 4.0)))


def chunks(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

class ConvBond(object):
    """docstring for ConvBond"""
    def __init__(self, bonds):
        super(ConvBond, self).__init__()
        # print bonds.types()
        # self.types = deepcopy(bonds.types())
        try:
            self.types = deepcopy(bonds.types())
            # print type(self.types)
        except IndexError, e:
            self.types = []
        self.indices = deepcopy(bonds.to_indices())

    def types(self):
        return self.types

    def to_indices(self):
        return self.indices


def loop_trajectory(sel1c,sel2c,indices1,indices2,config,suppl):
    # print(len(sel1c) , len(sel2c))
    # indices1 = suppl[0]
    # indices2 = suppl[1]
    cutoff, hbondcutoff, hbondcutangle = config
    # resname_array = comm.bcast(resname_array, root=0)
    # resid_array = comm.bcast(resid_array, root=0)
    # name_array = comm.bcast(name_array, root=0)
    type_array = suppl[0]
    bonds = suppl[1]
    # segids = comm.bcast(segids, root=0)
    # backbone = comm.bcast(backbone, root=0)
    heavyatoms = suppl[2]
    name_array = suppl[3]

    allRankContacts = []
    start = time.time()
    for s1, s2 in zip(sel1c, sel2c):
        frame = 0
        currentFrameContacts = []
        result = np.ndarray(shape=(len(s1), len(s2)), dtype=float)
        distarray = distances.distance_array(s1, s2, box=None, result=result)
        contacts = np.where(distarray <= cutoff)
        for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
            convindex1 = indices1[frame][idx1]  # idx1 converted to global atom indexing
            convindex2 = indices2[frame][idx2]  # idx2 converted to global atom indexing
            # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered, hydrogen bonds can still be detected!
            if re.match("H(.*)", name_array[convindex1]) or re.match("H(.*)", name_array[convindex2]):
                continue
            # distance between atom1 and atom2
            distance = distarray[idx1, idx2]
            weight = weight_function(distance)
            hydrogenBonds = []
            type1 = next((x.htype for x in heavyatoms if x.name == type_array[convindex1]), AtomHBondType.none)
            type2 = next((x.htype for x in heavyatoms if x.name == type_array[convindex2]), AtomHBondType.none)
            if type1 != AtomHBondType.none and type2 != AtomHBondType.none:
                if (type1 == AtomHBondType.both and type2 == AtomHBondType.both) or \
                        (type1 == AtomHBondType.acc and type2 == AtomHBondType.don) or \
                        (type1 == AtomHBondType.don and type2 == AtomHBondType.acc) or \
                        (type1 == AtomHBondType.both and type2 == AtomHBondType.acc) or \
                        (type1 == AtomHBondType.acc and type2 == AtomHBondType.both) or \
                        (type1 == AtomHBondType.don and type2 == AtomHBondType.both) or \
                        (type1 == AtomHBondType.both and type2 == AtomHBondType.don):
                    # print "hbond? %s - %s" % (type_array[convindex1], type_array[convindex2])
                    # search for hatom, check numbering in bond!!!!!!!!!!
                    b1 = bonds[convindex1]
                    b2 = bonds[convindex2]
                    # search for hydrogen atoms bound to atom 1
                    bondcount1 = 0
                    hydrogenAtomsBoundToAtom1 = []
                    for b in b1.types:
                        hydrogen = next((x for x in b if x.startswith("H")), 0)
                        # print b
                        if hydrogen != 0:
                            # print "h bond to atom1"
                            bondindices1 = b1.to_indices()[bondcount1]
                            hydrogenidx = next(
                                ((j + 1) for j in bondindices1 if type_array[j + 1].startswith("H")), -1)
                            if hydrogenidx != -1:
                                # print type_array[hydrogenidx]
                                hydrogenAtomsBoundToAtom1.append(hydrogenidx)
                        bondcount1 += 1
                    # search for hydrogen atoms bound to atom 2
                    bondcount2 = 0
                    hydrogenAtomsBoundToAtom2 = []
                    for b in b2.types:
                        hydrogen = next((x for x in b if x.startswith("H")), 0)
                        # print b
                        if hydrogen != 0:
                            # print "h bond to atom2"
                            bondindices2 = b2.to_indices()[bondcount2]
                            hydrogenidx = next(
                                ((k + 1) for k in bondindices2 if type_array[k + 1].startswith("H")), -1)
                            if hydrogenidx != -1:
                                # print type_array[hydrogenidx]
                                hydrogenAtomsBoundToAtom2.append(hydrogenidx)
                        bondcount2 += 1
                    # check hbond criteria for hydrogen atoms bound to first atom
                    for global_hatom in hydrogenAtomsBoundToAtom1:
                        conv_hatom = np.where(indices1[frame] == global_hatom)[0][0]
                        typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[convindex2]),
                                         AtomHBondType.none)
                        if typeHeavy == AtomHBondType.acc and (distarray[conv_hatom, idx2] <= hbondcutoff):
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
                                dist = distarray[conv_hatom, idx2]
                                new_hbond = HydrogenBond(convindex1, convindex2, global_hatom, dist, angle,
                                                         hbondcutoff,
                                                         hbondcutangle)
                                hydrogenBonds.append(new_hbond)
                            # print str(convindex1) + " " + str(convindex2)
                            # print "hbond found: %d,%d,%d"%(convindex1,global_hatom,convindex2)
                            # print angle
                    for global_hatom in hydrogenAtomsBoundToAtom2:
                        conv_hatom = np.where(indices2[frame] == global_hatom)[0][0]
                        typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[convindex1]),
                                         AtomHBondType.none)
                        if typeHeavy == AtomHBondType.acc and (distarray[idx1, conv_hatom] <= hbondcutoff):
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
                                dist = distarray[idx1, conv_hatom]
                                new_hbond = HydrogenBond(convindex2, convindex1, global_hatom, dist, angle,
                                                         hbondcutoff,
                                                         hbondcutangle)
                                hydrogenBonds.append(new_hbond)
            newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1), int(convindex2),
                                 hydrogenBonds)
            currentFrameContacts.append(newAtomContact)
        allRankContacts.append(currentFrameContacts)
        frame+=1
    return allRankContacts


def run_load_parallel(nproc, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
    # nproc = int(self.settingsView.coreBox.value())
    pool = LoggingPool(nproc)
    # manager = multiprocessing.Manager()
    # d=manager.list(trajArgs)

    heavyatomlines = []
    heavyatoms = []
    pars = open(os.path.dirname(os.path.abspath(__file__)) + '/testpar.prm', 'r')
    for line in pars:
        if re.match("MASS", line):
            heavyatomlines.append(line.rstrip())
    for atomline in heavyatomlines:
        # read new AtomType and its corresponding AtomHBondType from file
        atype = AtomType.parseParameterFileString(atomline)
        heavyatoms.append(atype)

    #load psf and dcd
    u = MDAnalysis.Universe(psf, dcd)
    # define selections according to sel1text and sel2text
    sel1 = u.select_atoms(sel1text)
    sel2 = u.select_atoms(sel2text)

    # write properties of all atoms to lists
    all_sel = u.select_atoms("all")
    backbone_sel = u.select_atoms("backbone")
    resname_array = []
    resid_array = []
    name_array = []
    type_array = []
    bonds = []
    segids = []
    backbone = []
    for atom in all_sel.atoms:
        resname_array.append(atom.resname)
        resid_array.append(atom.resid)
        name_array.append(atom.name)
        type_array.append(atom.type)
        bonds.append(ConvBond(atom.bonds))
        segids.append(atom.segid)
    for atom in backbone_sel:
        backbone.append(atom.index)

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

    contactResults = []
    # loop over trajectory
    totalFrameNumber = len(u.trajectory)
    start = time.time()
    sel1c = chunks(sel1coords, nproc)
    sel2c = chunks(sel2coords, nproc)
    sel1ind = chunks(indices1, nproc)
    sel2ind = chunks(indices2, nproc)
    print(len(sel1ind),len(sel2ind))
    # show trajectory information and selection information
    print("trajectory with %d frames loaded" % len(u.trajectory))

    print("Running on %d cores" % nproc)
    results = []
    rank = 0
    for c in zip(sel1c,sel2c,sel1ind,sel2ind):
        results.append( pool.apply_async( loop_trajectory, args=(c[0],c[1],c[2],c[3],[cutoff, hbondcutoff, hbondcutangle],[type_array,bonds,heavyatoms,name_array])) )
        rank +=1
    # TODO: might be important, but without, it's faster and until now, provides the same results
    pool.close()
    pool.join()
    stop = time.time()
    print("time: ", str(stop-start), rank)
    allContacts = []
    for res in results:
        rn = res.get()
        print(len(rn))
        allContacts.extend(rn)
    print("frames: ", len(allContacts))
    return [allContacts,resname_array,resid_array,name_array,type_array,segids,backbone,sel1text,sel2text]
