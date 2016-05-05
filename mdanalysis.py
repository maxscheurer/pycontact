#md analysis implementation for contact search
import MDAnalysis
import numpy as np
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.hbonds.hbond_analysis import *
import itertools
from timeit import default_timer as timer
import re
class AtomContact:
	def __init__(self, frame, distance, weight, idx1, idx2, hbondinfo, aType1, aType2, resid1, resid2, resname1, resname2, segid1, segid2):
		self.frame = int(frame)
		self.distance = float(distance)
		self.weight = float(weight)
		self.idx1 = int(idx1)
		self.idx2 = int(idx2)
		self.hbondinfo = hbondinfo
		self.aType1 = aType1
		self.aType2 = aType2
		self.resid1 = resid1
		self.resid2 = resid2
		self.resname1 = resname1
		self.resname2 = resname2
		self.segid1 = segid1
		self.segid2 = segid2
	def toString(self):
		print "frame: %d, dist: %f, weight: %f, idx1: %d, idx2: %d" % (self.frame,self.distance, self.weight, self.idx1, self.idx2)

class HydrogenBond:
	def __init__(self, donorIndex, acceptorIndex, hydrogenIndex, acceptorHydrogenDistance, angle, usedCutoffDist, usedCutoffAngle):
		self.donorIndex = donorIndex
		self.acceptorIndex = acceptorIndex
		self.hydrogenIndex = hydrogenIndex
		self.acceptorHydrogenDistance = acceptorHydrogenDistance
		self.angle = angle
		self.usedCutoffDist = usedCutoffDist
		self.usedCutoffAngle = usedCutoffAngle
	def toString(self):
		print "donor: %d, acceptor: %d, hydrogen: %d, dist: %f, angle: %f, usedCutoffDist: %f, usedCutoffAngle: %f" % (self.donorIndex, self.acceptorIndex, self.hydrogenIndex, self.acceptorHydrogenDistance, self.angle, self.usedCutoffDist, self.usedCutoffAngle)

def weight_function(value):
	return (1.0)/(1.0 + np.exp(5.0*(value-4.0)))

class AtomHBondType:
	don, acc, both, none = range(4)
	mapping = {"don":don, "acc":acc, "both":both, "none":none}

class AtomType:
	def __init__(self, name, comment, htype):
		self.name = name
		self.comment = comment
		self.htype = htype
	@staticmethod
	def parseParameterFileString(string):
		spl = string.split("!")
		name = spl[0].split()[2]
		comment = spl[1][1:]
		try:
			htype = spl[2]
		except:
			htype = "none"
		tp = AtomType(name,comment, AtomHBondType.mapping[htype])
		return tp

heavyatomlines = []
heavyatoms = []
pars = open('testpar.prm', 'r')
for line in pars:
    if re.match("MASS", line):
        heavyatomlines.append(line.rstrip())
for atomline in heavyatomlines:
	atype = AtomType.parseParameterFileString(atomline)
	heavyatoms.append(atype)

### config

cutoff = 5.0
hbondcutoff = 3.5
hbondcutangle = 120
sel1text = "segid RN11"
sel2text = "segid UBQ"

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
name_array = []
type_array = []
bonds = []
segids = []
for atom in all_sel.atoms:
	resname_array.append(atom.resname)
	resid_array.append(atom.resid)
	name_array.append(atom.name)
	type_array.append(atom.type)
	bonds.append(atom.bonds)
	segids.append(atom.segid)
print "trajectory with %d frames loaded" % len(u.trajectory)
print len(sel1.coordinates()),len(sel2.coordinates())
start = timer()
contactResults = []
#loop over trajectory
for ts in u.trajectory: 
	currentFrameContacts = []
	frame = ts.frame
	result = np.ndarray(shape=(len(sel1.coordinates()),len(sel2.coordinates())), dtype=float)
	distarray = distances.distance_array(sel1.coordinates(), sel2.coordinates(), box=None, result=result)
	contacts = np.where(distarray <= cutoff)
	for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
		convindex1 = indices1[idx1]
		convindex2 = indices2[idx2]
		#jump out of loop if hydrogen contacts are found
		if re.match("H(.*)", type_array[convindex1]) or re.match("H(.*)", type_array[convindex2]):
			continue  
		distance = distarray[idx1,idx2]
		weight = weight_function(distance)
		type1 = next((x.htype for x in heavyatoms if x.name == type_array[convindex1]), AtomHBondType.none)
		type2 = next((x.htype for x in heavyatoms if x.name == type_array[convindex2]), AtomHBondType.none)
		hydrogenBonds = []
		if type1 != AtomHBondType.none and type2 != AtomHBondType.none:
			if (type1 == AtomHBondType.both and type2 == AtomHBondType.both) or \
				(type1 == AtomHBondType.acc and type2 == AtomHBondType.don) or \
				(type1 == AtomHBondType.don and type2 == AtomHBondType.acc) or \
				(type1 == AtomHBondType.both and type2 == AtomHBondType.acc) or \
				(type1 == AtomHBondType.acc and type2 == AtomHBondType.both) or \
				(type1 == AtomHBondType.don and type2 == AtomHBondType.both) or \
				(type1 == AtomHBondType.both and type2 == AtomHBondType.don):
					# print "hbond? %s - %s" % (type_array[convindex1], type_array[convindex2])
					#search for hatom, check numbering in bond!!!!!!!!!!
					b1 = bonds[convindex1]
					b2 = bonds[convindex2]
					# search for hydrogen atoms bound to atom 1
					bondcount1 = 0
					hydrogenAtomsBoundToAtom1 = []
					for b in b1.types():
						hydrogen = next((x for x in b if x.startswith("H")), 0)
						# print b
						if hydrogen != 0:
							# print "h bond to atom1"
							bondindices1 = b1.to_indices()[bondcount1]
							hydrogenidx = next(((j+1) for j in bondindices1 if type_array[j+1].startswith("H")), -1)
							if hydrogenidx != -1:
								# print type_array[hydrogenidx]
								hydrogenAtomsBoundToAtom1.append(hydrogenidx)
						bondcount1 += 1
					# search for hydrogen atoms bound to atom 2
					bondcount2 = 0
					hydrogenAtomsBoundToAtom2 = []
					for b in b2.types():
						hydrogen = next((x for x in b if x.startswith("H")), 0)
						# print b
						if hydrogen != 0:
							# print "h bond to atom2"
							bondindices2 = b2.to_indices()[bondcount2]
							hydrogenidx = next(((k+1) for k in bondindices2 if type_array[k+1].startswith("H")), -1)
							if hydrogenidx != -1:
								# print type_array[hydrogenidx]
								hydrogenAtomsBoundToAtom2.append(hydrogenidx)
						bondcount2 += 1
					# check hbond criteria for hydrogen atoms bound to first atom
					for global_hatom in hydrogenAtomsBoundToAtom1:
						conv_hatom = indices1.index(global_hatom)
						scan = np.where(distarray[conv_hatom,:] <= hbondcutoff)
						for idxheavy in scan[0]:
							heavy2_converted = indices2[idxheavy]
							typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[heavy2_converted]), AtomHBondType.none)
							if typeHeavy == AtomHBondType.acc:
								donorPosition = sel1.coordinates()[idx1]
								hydrogenPosition = sel1.coordinates()[conv_hatom]
								acceptorPosition = sel2.coordinates()[idxheavy]
								v1 = hydrogenPosition - acceptorPosition
								v2 = hydrogenPosition - donorPosition
								v1norm = np.linalg.norm(v1)
								v2norm = np.linalg.norm(v2)
								dot = np.dot(v1,v2)
								angle = np.degrees(np.arccos(dot/(v1norm*v2norm)))
								if angle >= hbondcutangle:
									dist = distarray[conv_hatom,idxheavy]
									new_hbond = HydrogenBond(convindex1,heavy2_converted,global_hatom, dist, angle, hbondcutoff, hbondcutangle)
									hydrogenBonds.append(new_hbond)
									# print "hbond found: %d,%d,%d"%(convindex1,global_hatom,heavy2_converted)
									# print angle
					for global_hatom in hydrogenAtomsBoundToAtom2:
						conv_hatom = indices2.index(global_hatom)
						scan = np.where(distarray[:,conv_hatom] <= hbondcutoff)
						for idxheavy in scan[0]:
							heavy1_converted = indices1[idxheavy]
							typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[heavy1_converted]), AtomHBondType.none)
							if typeHeavy == AtomHBondType.acc:
								donorPosition = sel2.coordinates()[idx2]
								hydrogenPosition = sel2.coordinates()[conv_hatom]
								acceptorPosition = sel1.coordinates()[idxheavy]
								v1 = hydrogenPosition - acceptorPosition
								v2 = hydrogenPosition - donorPosition
								v1norm = np.linalg.norm(v1)
								v2norm = np.linalg.norm(v2)
								dot = np.dot(v1,v2)
								angle = np.degrees(np.arccos(dot/(v1norm*v2norm)))
								if angle >= hbondcutangle:
									dist = distarray[idxheavy,conv_hatom]
									new_hbond = HydrogenBond(convindex2,heavy1_converted,global_hatom, dist, angle, hbondcutoff, hbondcutangle)
									hydrogenBonds.append(new_hbond)
									# print "hbond found: %d,%d,%d"%(convindex2,global_hatom,heavy1_converted)
									# print angle
		#finalize
		newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1), int(convindex2),hydrogenBonds, type_array[convindex1], type_array[convindex2], resid_array[convindex1], resid_array[convindex2], resname_array[convindex1], resname_array[convindex2], segids[convindex1], segids[convindex2])
		currentFrameContacts.append(newAtomContact)
	contactResults.append(currentFrameContacts)
#draft & sketch playground

#prototype writer for vmd visualization
# f = open('showHBondsInVMD.tcl', 'w')
# f.write('mol new %s \n'%psf)
# f.write('mol addfile %s \n'%dcd)
# f.write('mol delrep 0 top \n')
# f.write('mol representation NewCartoon \n')
# f.write('mol Color ColorID 3 \n')
# f.write('mol selection {all} \n')
# f.write('mol addrep top \n')
# for frame in contactResults:
# 	for contact in frame:
# 		for hbond in contact.hbondinfo:
# 			# hbond.toString()
# 			f.write('mol representation VDW \n')
# 			f.write('mol Color Name \n')
# 			f.write('mol selection {index %d %d %d} \n'%(hbond.donorIndex, hbond.acceptorIndex, hbond.hydrogenIndex))
			# f.write('mol addrep top \n')
# f.close()

#frame analysis
contact_accumulated = []
allkeys = []
for frame in contactResults:
	currentFrameAcc = {}
	for cont in frame:
		key = '%s %s %s %s' % (cont.resid1, cont.resname1, cont.resid2, cont.resname2)
		if key in currentFrameAcc:
			currentFrameAcc[key]+=cont.weight
		else:
			currentFrameAcc[key]=cont.weight
		if not key in allkeys:
			allkeys.append(key)
	contact_accumulated.append(currentFrameAcc)

for key in allkeys:
	for frame_dict in contact_accumulated:
		if not key in frame_dict:
			frame_dict[key] = 0
		# print frame_dict
# print analysis time and quit
stop = timer()
print (stop - start)
quit()
