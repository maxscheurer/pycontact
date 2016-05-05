#md analysis implementation for contact search
import MDAnalysis
import numpy as np
from MDAnalysis.analysis import distances
import itertools
from timeit import default_timer as timer

class AtomContact:
	def __init__(self, frame, distance, weight, idx1, idx2):
		self.frame = int(frame)
		self.distance = float(distance)
		self.weight = float(weight)
		self.idx1 = int(idx1)
		self.idx2 = int(idx2)
	def toString(self):
		print "frame: %d, dist: %f, weight: %f, idx1: %d, idx2: %d" % (self.frame,self.distance, self.weight, self.idx1, self.idx2)

cutoff = 5.0
sel1text = "segid RN11 and not name H*"
sel2text = "segid UBQ and not name H*"

psf = "rpn11_ubq_interface-ionized.psf" 
dcd = "short.dcd"

### main tool
u = MDAnalysis.Universe(psf,dcd) 
sel1 = u.select_atoms(sel1text)
sel2 = u.select_atoms(sel2text)
indices1 = []
for at in sel1.atoms:
	indices1.append(at.index)
indices2 = []
for at in sel2.atoms:
	indices2.append(at.index)
all_sel = u.select_atoms("all")
resname_array = []
resid_array = []
for atom in all_sel.atoms:
	resname_array.append(atom.resname)
	resid_array.append(atom.resid)
print "trajectory with %d frames loaded" % len(u.trajectory)
print len(sel1.coordinates()),len(sel2.coordinates())
start = timer()
contactResults = []
for ts in u.trajectory: 
	frame = ts.frame
	result = np.ndarray(shape=(len(sel1.coordinates()),len(sel2.coordinates())), dtype=float)
	distarray = distances.distance_array(sel1.coordinates(), sel2.coordinates(), box=None, result=result)
	contacts = np.where(distarray <= cutoff)
	for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
		convindex1 = indices1[idx1]
		convindex2 = indices2[idx2]
		distance = distarray[idx1,idx2]
		weight = 0
		newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1), int(convindex2))
		contactResults.append(newAtomContact)
		# print resname_array[convindex1] + "-" + str(resid_array[convindex1]) + "  " + resname_array[convindex2] + "-" + str(resid_array[convindex2]) + " " + str(distarray[idx1,idx2])
stop = timer()
for contact in contactResults:
	contact.toString()
print (stop - start)
quit()