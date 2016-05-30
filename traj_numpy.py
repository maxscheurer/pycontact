import MDAnalysis
import time
from mpi4py import MPI
from MDAnalysis.analysis import distances
import numpy as np
from mdanalysis import *
import itertools
import re, sys
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

cutoff = 5

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


if rank == 0:
	# load psf and dcd file in memory
	u = MDAnalysis.Universe("rpn11_ubq_interface-ionized.psf", "short.dcd")
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
	    bonds.append(atom.bonds)
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
		print ts.frame
		sel1coords.append(sel1.coordinates())
		sel2coords.append(sel2.coordinates())
	sel1c = chunks(sel1coords,size)
	sel2c = chunks(sel2coords,size)
else:
	print rank
	sel1c = None
	sel2c = None
	indices1 = None
	indices2 = None

indices1 = comm.bcast(indices1,root=0)
indices2 = comm.bcast(indices2,root=0)
sel1c = comm.scatter(sel1c,root=0)
sel2c = comm.scatter(sel2c,root=0)

for s1,s2 in zip(sel1c,sel2c):
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
        newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1),int(convindex2),hydrogenBonds)
        currentFrameContacts.append(newAtomContact)


# contactResults = []
# # loop over trajectory
# totalFrameNumber = len(u.trajectory)
# stop = time.time()
# print stop - start