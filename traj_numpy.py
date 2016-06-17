import MDAnalysis
import time
from mpi4py import MPI
from MDAnalysis.analysis import distances
import numpy as np
from mdanalysis import *
import itertools
import re, sys
from copy import deepcopy
import time



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


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

cutoff = 5
hbondcutoff = 2.5
hbondcutangle = 120

if rank == 0:
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
    # load psf and dcd file in memory
    
    u = MDAnalysis.Universe("rpn11_ubq_interface-ionized.psf", "short.dcd")
    # u = MDAnalysis.Universe("/Users/maximilianscheurer/Dropbox/TCBG/ba/data/prot_ubp6/yeast_proteasome_ubp6_corr.psf", "/Users/maximilianscheurer/Dropbox/TCBG/ba/data/prot_ubp6/10ns_mdff_prot_wubp6.dcd")
    # define selections according to sel1text and sel2text
    sel1 = u.select_atoms("segid RN11")
    sel2 = u.select_atoms("segid UBQ")
    # write atomindices for each selection to list
    indices1 = []
    for at in sel1.atoms:
        indices1.append(at.index)
    indices2 = []
    for at in sel2.atoms:
        indices2.append(at.index)
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
    # show trajectory information and selection information
    print "trajectory with %d frames loaded" % len(u.trajectory)
    print len(sel1.coordinates()), len(sel2.coordinates())
    sel1coords = []
    sel2coords = []
    start = time.time()
    for ts in u.trajectory:
        sel1coords.append(sel1.coordinates())
        sel2coords.append(sel2.coordinates())
    sel1c = chunks(sel1coords, size)
    sel2c = chunks(sel2coords, size)
else:
    print rank
    sel1c = None
    sel2c = None
    indices1 = None
    indices2 = None
    resname_array = None
    resid_array = None
    name_array = None
    type_array = None
    bonds = None
    segids = None
    backbone = None
    heavyatoms = None

indices1 = comm.bcast(indices1, root=0)
indices2 = comm.bcast(indices2, root=0)

# resname_array = comm.bcast(resname_array, root=0)
# resid_array = comm.bcast(resid_array, root=0)
# name_array = comm.bcast(name_array, root=0)
type_array = comm.bcast(type_array, root=0)
bonds = comm.bcast(bonds, root=0)
# segids = comm.bcast(segids, root=0)
# backbone = comm.bcast(backbone, root=0)
heavyatoms = comm.bcast(heavyatoms, root=0)

sel1c = comm.scatter(sel1c, root=0)
sel2c = comm.scatter(sel2c, root=0)

allRankContacts = []
start = time.time()
for s1, s2 in zip(sel1c, sel2c):
    frame = 0
    currentFrameContacts = []
    result = np.ndarray(shape=(len(s1), len(s2)), dtype=float)
    distarray = distances.distance_array(s1, s2, box=None, result=result)
    contacts = np.where(distarray <= cutoff)
    for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
        convindex1 = indices1[idx1]  # idx1 converted to global atom indexing
        convindex2 = indices2[idx2]  # idx2 converted to global atom indexing
        # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered, hydrogen bonds can still be detected!
        # if re.match("H(.*)", self.name_array[convindex1]) or re.match("H(.*)", self.name_array[convindex2]):
        #     continue
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
                    conv_hatom = indices1.index(global_hatom)
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
                    conv_hatom = indices2.index(global_hatom)
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

stop = time.time()
print rank,stop-start
allChunks = comm.gather(allRankContacts,root=0)
if rank == 0:
    allContacts = []
    for chunk in allChunks:
        for el in chunk:
            allContacts.append(el)
    print len(allContacts)
# contactResults = []
# # loop over trajectory
# totalFrameNumber = len(u.trajectory)
# stop = time.time()
# print stop - start
