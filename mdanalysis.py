#md analysis implementation for contact search
import MDAnalysis
import numpy as np
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.hbonds.hbond_analysis import *
import itertools
from timeit import default_timer as timer
import re
class AtomContact:
	def __init__(self, frame, distance, weight, idx1, idx2, hbondinfo):
		self.frame = int(frame)
		self.distance = float(distance)
		self.weight = float(weight)
		self.idx1 = int(idx1)
		self.idx2 = int(idx2)
		self.hbondinfo = hbondinfo
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

class AccumulatedContact(object):
	"""docstring for AccumulatedContact"""
	def __init__(self, key1,key2):
		super(AccumulatedContact, self).__init__()
		self.scoreArray = []
		self.contributingAtoms = []
		self.key1 = key1
		self.key2 = key2
	def addScore(self,newScore):
		self.scoreArray.append(newScore)
	def addContributingAtoms(self, contAtoms):
		self.contributingAtoms.append(contAtoms)


class TempContactAccumulate(object):
	"""docstring for TempContactAccumulate"""
	def __init__(self, key1, key2):
		super(TempContactAccumulate, self).__init__()
		self.fscore = 0
		self.contributingAtomContacts = []
		self.key1 = key1
		self.key2 = key2
		
class AccumulationMapIndex():
	index, atype, name, resid, resname, segid = range(6)
	mapping = ["i.", "t.", "nm.", "r.", "rn.", "s."]
		
def makeKeyArraysFromMaps(map1,map2,contact):
	idx1 = contact.idx1
	idx2 = contact.idx2
	counter = 0
	keys1 = []
	for val in map1:
		if val == 1:	
			if counter == AccumulationMapIndex.index:
				keys1.append(idx1)
			elif counter == AccumulationMapIndex.atype:
				keys1.append(type_array[idx1])
			elif counter == AccumulationMapIndex.name:
				keys1.append(name_array[idx1])
			elif counter == AccumulationMapIndex.resid:
				keys1.append(resid_array[idx1])
			elif counter == AccumulationMapIndex.resname:
				keys1.append(resname_array[idx1])
			elif counter == AccumulationMapIndex.segid:
				keys1.append(segids[idx1])
		else:
			keys1.append("none")		
		counter += 1
	counter = 0
	keys2 = []
	for val in map2:
		if val == 1:	
			if counter == AccumulationMapIndex.index:
				keys2.append(idx2)
			elif counter == AccumulationMapIndex.atype:
				keys2.append(type_array[idx2])
			elif counter == AccumulationMapIndex.name:
				keys2.append(name_array[idx2])
			elif counter == AccumulationMapIndex.resid:
				keys2.append(resid_array[idx2])
			elif counter == AccumulationMapIndex.resname:
				keys2.append(resname_array[idx2])
			elif counter == AccumulationMapIndex.segid:
				keys2.append(segids[idx2])
		else:
			keys2.append("none")		
		counter += 1
	return [keys1,keys2]

def makeKeyArraysFromKey(key):
	keystring1, keystring2 = key.split("-")
	mapping = AccumulationMapIndex.mapping
	maximal = len(mapping)
	key1 = []
	for i in range(0, maximal):
		current = mapping[i]
		if current not in keystring1:
			key1.append("none")
			continue
		if i == (maximal -1):
			key1.append(keystring1[keystring1.index(current) + len(current):])
			break
		nextCurrent = mapping[i+1]
		if nextCurrent not in keystring1:
			nxt = ""
			for k in mapping[i+1:]:
				if k in keystring1[keystring1.index(current) + len(current):]:
					nxt = k
					break
			if nxt != "":
				key1.append(keystring1[keystring1.index(current) + len(current):keystring1.index(nxt)])
			else:
				key1.append(keystring1[keystring1.index(current) + len(current):])
			continue
		else:
			currentValue = find_between(keystring1,current,nextCurrent)
			if currentValue == "":
				key1.append("none")
			else:
				key1.append(currentValue)

	key2 = []
	for i in range(0, maximal):
		current = mapping[i]
		if current not in keystring2:
			key2.append("none")
			continue
		if i == (maximal -1):
			key2.append(keystring2[keystring2.index(current) + len(current):])
			break
		nextCurrent = mapping[i+1]
		if nextCurrent not in keystring2:
			nxt = ""
			for k in mapping[i+1:]:
				if k in keystring2[keystring2.index(current) + len(current):]:
					nxt = k
					break
			if nxt != "":
				key2.append(keystring2[keystring2.index(current) + len(current):keystring2.index(nxt)])
			else:
				key2.append(keystring2[keystring2.index(current) + len(current):])
			continue
		else:
			currentValue = find_between(keystring2,current,nextCurrent)
			if currentValue == "":
				key2.append("none")
			else:
				key2.append(currentValue)
	return [key1,key2]


def makeKeyFromKeyArrays(key1,key2):
	key = ""
	itemcounter = 0
	for item in key1:
		if item != "none":
			key += AccumulationMapIndex.mapping[itemcounter] + str(item)
		itemcounter += 1
	key += "-"
	itemcounter = 0
	for item in key2:
		if item != "none":
			key += AccumulationMapIndex.mapping[itemcounter] + str(item)
		itemcounter += 1
	return key

heavyatomlines = []
heavyatoms = []
pars = open('testpar.prm', 'r')
for line in pars:
    if re.match("MASS", line):
        heavyatomlines.append(line.rstrip())
for atomline in heavyatomlines:
	atype = AtomType.parseParameterFileString(atomline)
	heavyatoms.append(atype)

def find_between(s, first, last):
    try:
        start = s.index(first) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

### config

cutoff = 5.0
hbondcutoff = 3.5
hbondcutangle = 120
sel1text = "segid RN11"
sel2text = "segid UBQ"

psf = "rpn11_ubq_interface-ionized.psf" 
dcd = "short.dcd"
# dcd = "rpn11_ubq_50ns.dcd"
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
	print frame
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
						typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[convindex2]), AtomHBondType.none)
						if typeHeavy == AtomHBondType.acc and (distarray[conv_hatom,idx2] <= hbondcutoff):
							donorPosition = sel1.coordinates()[idx1]
							hydrogenPosition = sel1.coordinates()[conv_hatom]
							acceptorPosition = sel2.coordinates()[idx2]
							v1 = hydrogenPosition - acceptorPosition
							v2 = hydrogenPosition - donorPosition
							v1norm = np.linalg.norm(v1)
							v2norm = np.linalg.norm(v2)
							dot = np.dot(v1,v2)
							angle = np.degrees(np.arccos(dot/(v1norm*v2norm)))
							if angle >= hbondcutangle:
								dist = distarray[conv_hatom,idx2]
								new_hbond = HydrogenBond(convindex1,convindex2,global_hatom, dist, angle, hbondcutoff, hbondcutangle)
								hydrogenBonds.append(new_hbond)
								# print str(convindex1) + " " + str(convindex2)
								# print "hbond found: %d,%d,%d"%(convindex1,global_hatom,convindex2)
								# print angle
					for global_hatom in hydrogenAtomsBoundToAtom2:
						conv_hatom = indices2.index(global_hatom)
						typeHeavy = next((x.htype for x in heavyatoms if x.name == type_array[convindex1]), AtomHBondType.none)
						if typeHeavy == AtomHBondType.acc and (distarray[idx1,conv_hatom] <= hbondcutoff):
							donorPosition = sel2.coordinates()[idx2]
							hydrogenPosition = sel2.coordinates()[conv_hatom]
							acceptorPosition = sel1.coordinates()[idx1]
							v1 = hydrogenPosition - acceptorPosition
							v2 = hydrogenPosition - donorPosition
							v1norm = np.linalg.norm(v1)
							v2norm = np.linalg.norm(v2)
							dot = np.dot(v1,v2)
							angle = np.degrees(np.arccos(dot/(v1norm*v2norm)))
							if angle >= hbondcutangle:
								dist = distarray[idx1,conv_hatom]
								new_hbond = HydrogenBond(convindex2,convindex1,global_hatom, dist, angle, hbondcutoff, hbondcutangle)
								hydrogenBonds.append(new_hbond)
								# print str(convindex1) + " " + str(convindex2)
								# print "hbond found: %d,%d,%d"%(convindex2,global_hatom,convindex1)
								# print angle
		#finalize
		# , type_array[convindex1], type_array[convindex2], resid_array[convindex1], resid_array[convindex2], resname_array[convindex1], resname_array[convindex2], segids[convindex1], segids[convindex2]
		newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1), int(convindex2),hydrogenBonds)
		currentFrameContacts.append(newAtomContact)
	contactResults.append(currentFrameContacts)

map1 = [0,0,0,1,1,0]
map2 = [0,0,0,1,1,0]
frame_contacts_accumulated = []
allkeys = []
for frame in contactResults:
	currentFrameAcc = {}
	for cont in frame:
		key1, key2 = makeKeyArraysFromMaps(map1,map2,cont)
		key = makeKeyFromKeyArrays(key1, key2)
		if key in currentFrameAcc:
			currentFrameAcc[key].fscore+=cont.weight
			currentFrameAcc[key].contributingAtomContacts.append(cont)
		else:
			currentFrameAcc[key] = TempContactAccumulate(key1, key2)
			currentFrameAcc[key].fscore+=cont.weight
			currentFrameAcc[key].contributingAtomContacts.append(cont)
		if not key in allkeys:
			allkeys.append(key)
	frame_contacts_accumulated.append(currentFrameAcc)

accumulatedContacts = {}
for key in allkeys:
	accumulatedContacts[key] = []
	for frame_dict in frame_contacts_accumulated:
		if not key in frame_dict:
			key1, key2 = makeKeyArraysFromKey(key)
			emptyCont = TempContactAccumulate(key1,key2)
			emptyCont.fscore = 0
			frame_dict[key] = emptyCont
		accumulatedContacts[key].append(frame_dict[key])

finalAccumulatedContacts = []
for key in accumulatedContacts:
	key1, key2 = makeKeyArraysFromKey(key)
	acc = AccumulatedContact(key1,key2)
	for tempContact in accumulatedContacts[key]:
		acc.addScore(tempContact.fscore)
		acc.addContributingAtoms(tempContact.contributingAtomContacts)
	finalAccumulatedContacts.append(acc)

# print analysis time and quit
stop = timer()
print (stop - start)
quit()