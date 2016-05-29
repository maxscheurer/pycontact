'''
    Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
'''
# md analysis implementation for contact search
# May 2016
# Author: Maximilian Scheurer, mscheurer@ks.uiuc.edu
import MDAnalysis
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.hbonds.hbond_analysis import *
import itertools
import re, sys
from biochemistry import *
from read_db import *
from copy import deepcopy
import time

MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

# type of Atom concerning it hbond behavior
class AtomHBondType:
    don, acc, both, none = range(4)
    mapping = {"don": don, "acc": acc, "both": both, "none": none}


# AtomType represents MASS entry for atoms in CHARMM topology and parameter files
class AtomType:
    def __init__(self, name, comment, htype):
        self.name = name  # name = atomtype in CHARMM file
        self.comment = comment  # properties/infos according to CHARMM file
        self.htype = htype  # AtomHBondType of AtomType, parsed from extra comment in CHARMM file

    @staticmethod
    def parseParameterFileString(string):
        # read charmm parameter/topology file to determine AtomTypes and their AtomHBondType
        spl = string.split("!")
        name = spl[0].split()[2]
        comment = spl[1][1:]
        try:
            htype = spl[2]
        except:
            htype = "none"
        tp = AtomType(name, comment, AtomHBondType.mapping[htype])
        return tp


## contains all information of a contact, mapped by key1/key2
# Main class for future implementation in PyContact
# will replace the Contact class in PyContact
class AccumulatedContact(object):
    """docstring for AccumulatedContact"""

    def __init__(self, key1, key2):
        super(AccumulatedContact, self).__init__()
        self.scoreArray = []  # score for every frame, summed by settings given in key1,key2
        self.contributingAtoms = []  # list of list of AtomContacts, contributing to the contact
        self.key1 = key1  # list of properties of sel1, cf. AccumulationMapIndex
        self.key2 = key2  # list of properties of sel2, cf. AccumulationMapIndex
        # TODO:  human readable implementation of key/title
        # self.title = makeKeyFromKeyArrays(key1, key2)
        self.title = self.human_readable_title()
        self.bb1 = 0
        self.sc1 = 0
        self.bb2 = 0
        self.sc2 = 0

    def human_readable_title(self):
        total = []
        for key in [self.key1, self.key2]:
            titleDict = {}
            counter = 0
            for item in key:
                titleDict[AccumulationMapIndex.mapping[counter]] = (item if item != "none" else "")
                counter += 1
            residueString = "%s%s" % (titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.resname]],
                                      str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.resid]]))
            atomIndexString = ("%s %s" % (AccumulationMapIndex.mapping[AccumulationMapIndex.index],
                                          str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.index]])) if
                               titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.index]] != "" else "")
            atomNameString = ("%s %s" % (AccumulationMapIndex.mapping[AccumulationMapIndex.name],
                                         str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.name]])) if
                              titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.name]] != "" else "")
            segnameString = ("%s %s" % (AccumulationMapIndex.mapping[AccumulationMapIndex.segid],
                                        str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.segid]])) if
                             titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.segid]] != "" else "")
            list = [residueString, atomIndexString, atomNameString, segnameString]
            finishedList = []
            for string in list:
                if string != "":
                    finishedList.append(string)
            finishedString = " , ".join(finishedList)
            total.append(finishedString)
        return " - ".join(total)

    def addScore(self, newScore):
        # append a score to the scoreArray, e.g. when a new frame score is added
        self.scoreArray.append(newScore)

    def determineBackboneSidechainType(self):
        if self.bb1 > self.sc1:
            self.atom1contactsBy = BackboneSidechainType.contactsBb
        else:
            self.atom1contactsBy = BackboneSidechainType.contactsSc

        if self.bb2 > self.sc2:
            self.atom2contactsBy = BackboneSidechainType.contactsBb
        else:
            self.atom2contactsBy = BackboneSidechainType.contactsSc

        if self.atom1contactsBy == BackboneSidechainType.contactsBb and self.atom2contactsBy == BackboneSidechainType.contactsBb:
            self.backboneSideChainType = BackboneSidechainContactType.bb_only
        elif self.atom1contactsBy == BackboneSidechainType.contactsSc and self.atom2contactsBy == BackboneSidechainType.contactsSc:
            self.backboneSideChainType = BackboneSidechainContactType.sc_only
        else:
            self.backboneSideChainType = BackboneSidechainContactType.both
        return self.backboneSideChainType

    def addContributingAtoms(self, contAtoms):
        # append a list of contributing atom to the contributingAtoms list, e.g. when a new frame is added
        self.contributingAtoms.append(contAtoms)  ## used for temporary accumulation of contacts in data analysis

    def setScores(self):
    	self.mean_score()
    	self.median_score()

    def hbond_percentage(self):
        self.hbondFramesScan()
        fnumber = len(self.scoreArray)
        counter = 0
        for element in self.hbondFrames:
            if element > 0:
                counter += 1
        return float(counter)/float(fnumber) * 100

    def total_time(self, ns_per_frame, threshold):
        time = 0
        for score in self.scoreArray:
            if score > threshold:
                time += ns_per_frame
        self.ttime = time
        return self.ttime

    def mean_life_time(self, ns_per_frame, threshold):
        self.meanLifeTime = np.mean(self.life_time(ns_per_frame, threshold))
        return self.meanLifeTime

    def median_life_time(self, ns_per_frame, threshold):
        self.medianLifeTime = np.median(self.life_time(ns_per_frame, threshold))
        return self.medianLifeTime

    def mean_score(self):
        mean = 0
        for score in self.scoreArray:
            mean += score
        mean = mean / len(self.scoreArray)
        self.meanScore = mean
        return mean

    def median_score(self):
        med = np.median(self.scoreArray)
        self.medianScore = med
        return med

    def life_time(self, ns_per_frame, threshold):
        lifeTimes = []
        contactActive = False
        contactTime = 0
        i = 0
        for score in self.scoreArray:
            if contactActive == False and score > threshold:
                contactActive = True
                contactTime += ns_per_frame
            elif contactActive == True and score > threshold:
                contactTime += ns_per_frame
            elif contactActive == True and score <= threshold:
                contactActive = False
                lifeTimes.append(contactTime)
                contactTime = 0
            if i == (len(self.scoreArray) - 1):
                lifeTimes.append(contactTime)
            i += 1
        return lifeTimes

    def hbondFramesScan(self):
        self.hbondFrames = []
        for frameList in self.contributingAtoms:
            currentFrame = 0
            for contAtoms in frameList:
                length = len(contAtoms.hbondinfo)
                if length > 0:
                    currentFrame += length
            self.hbondFrames.append(currentFrame)
        return self.hbondFrames

    def setContactType(self):
        self.contactType = self.determine_ctype()

    def determine_ctype(self):
        #only works if both maps contain resname
        r1 = self.key1[AccumulationMapIndex.resname].lower()
        r2 = self.key2[AccumulationMapIndex.resname].lower()
        if r1 == "none" or r2 == "none":
            return ContactType.other
        self.determineBackboneSidechainType()
        try:
            sc1 = str(read_residue_db("scpolarity", "name", r1)[0]["scpolarity"])
            scpol1 = SideChainPolarity.mapping[sc1]
        except:
            scpol1 = SideChainPolarity.other

        try:
            sc2 = str(read_residue_db("scpolarity", "name", r2)[0]["scpolarity"])
            scpol2 = SideChainPolarity.mapping[sc2]
        except:
            scpol2 = SideChainPolarity.other

        # hydrogen bonds: donor, acceptor, both

        try:
            hb1 = str(read_residue_db("hbondtype", "name", r1)[0]["hbondtype"])
            hbond1 = HBondType.mapping[hb1]
        except:
            hbond1 = HBondType.other

        try:
            hb2 = str(read_residue_db("hbondtype", "name", r2)[0]["hbondtype"])
            hbond2 = HBondType.mapping[hb2]
        except:
            hbond2 = HBondType.other

        # check if both residues contact by backbone
        if self.atom1contactsBy == BackboneSidechainType.contactsBb and self.atom2contactsBy == BackboneSidechainType.contactsBb:
            return ContactType.hbond
        # check if contact by sidechain and backbone
        if (self.atom1contactsBy == BackboneSidechainType.contactsBb and self.atom2contactsBy == BackboneSidechainType.contactsSc):
            # sc contact by residue B
            if hbond2 != HBondType.none:
                return ContactType.hbond
            else:
                return ContactType.other
        elif (self.atom1contactsBy == BackboneSidechainType.contactsSc and self.atom2contactsBy == BackboneSidechainType.contactsBb):
            # sc contact by residue A
            if hbond1 != HBondType.none:
                return ContactType.hbond
            else:
                return ContactType.other
        if self.atom1contactsBy == BackboneSidechainType.contactsSc and self.atom2contactsBy == BackboneSidechainType.contactsSc:
            # check for saltbridge
            if (scpol1 == SideChainPolarity.positive and scpol2 == SideChainPolarity.negative) or \
                    (scpol2 == SideChainPolarity.positive and scpol1 == SideChainPolarity.negative):
                return ContactType.saltbr
            # check for hydrophobic contact
            if scpol1 == SideChainPolarity.nonpolar and scpol2 == SideChainPolarity.nonpolar:
                if hbond1 == HBondType.none or hbond2 == HBondType.none:
                    return ContactType.hydrophobic
            # final hbond scan
            if (hbond1 == HBondType.donor and hbond2 == HBondType.acceptor) or \
                    (hbond1 == HBondType.acceptor and hbond2 == HBondType.donor) or \
                    (hbond1 == HBondType.both and hbond2 == HBondType.both) or \
                    (hbond1 == HBondType.donor and hbond2 == HBondType.both) or \
                    (hbond1 == HBondType.both and hbond2 == HBondType.donor) or \
                    (hbond1 == HBondType.acceptor and hbond2 == HBondType.both) or \
                    (hbond1 == HBondType.both and hbond2 == HBondType.acceptor):
                return ContactType.hbond

            return ContactType.other


# stores the frame's score as well as the key
# many TempContactAccumulated objects are later converted to AccumulatedContact
class TempContactAccumulate(object):
    """docstring for TempContactAccumulate"""

    def __init__(self, key1, key2):
        super(TempContactAccumulate, self).__init__()
        self.fscore = 0  # score of current frame
        self.contributingAtomContacts = []  # contrib. atoms, later appended to AccumulatedContact's contributingAtoms list
        self.key1 = key1
        self.key2 = key2
        self.bb1score = 0
        self.bb2score = 0
        self.sc1score = 0
        self.sc2score = 0


## enum and mapping for atom properties
# used to dynamically define keys and have a bijective nomenclature for all properties
class AccumulationMapIndex():
    index, atype, name, resid, resname, segid = range(6)
    mapping = ["i.", "t.", "nm.", "r.", "rn.", "s."]
    vmdsel = ["index", "type", "name", "resid", "resname", "segname"]


# contains infos about an atom-atom contact
class AtomContact:
    def __init__(self, frame, distance, weight, idx1, idx2, hbondinfo):
        self.frame = int(frame)  # frame the contact occured in
        self.distance = float(distance)  # distance between atom1 and atom2
        self.weight = float(weight)  # weighted distance according to applied weight function
        self.idx1 = int(idx1)  # global!!! index of atom1
        self.idx2 = int(idx2)  # global!!! index of atom2
        self.hbondinfo = hbondinfo  # list of HydrogenBonds for frame, empty list if no hbonds occured

    def toString(self):
        # print details about contact to console
        print "frame: %d, dist: %f, weight: %f, idx1: %d, idx2: %d" % (
            self.frame, self.distance, self.weight, self.idx1, self.idx2)


# contains infos about a hydrogenbond corresponding to a contact between two heavyatoms
class HydrogenBond:
    def __init__(self, donorIndex, acceptorIndex, hydrogenIndex, acceptorHydrogenDistance, angle, usedCutoffDist,
                 usedCutoffAngle):
        self.donorIndex = donorIndex  # global! index of hbond donor: D-h ... a
        self.acceptorIndex = acceptorIndex  # global! index of acceptor: d-h ... A
        self.hydrogenIndex = hydrogenIndex  # global! index of hydrogen atom: d-H ... a
        self.acceptorHydrogenDistance = acceptorHydrogenDistance  # distance of H and A: H-(acceptorHydrogenDistance)-A
        self.angle = angle  # angle between D-H-A
        self.usedCutoffDist = usedCutoffDist  # obvious
        self.usedCutoffAngle = usedCutoffAngle  # obvious

    def toString(self):
        # print details about hbond to console
        print "donor: %d, acceptor: %d, hydrogen: %d, dist: %f, angle: %f, usedCutoffDist: %f, usedCutoffAngle: %f" % (
            self.donorIndex, self.acceptorIndex, self.hydrogenIndex, self.acceptorHydrogenDistance, self.angle,
            self.usedCutoffDist, self.usedCutoffAngle)


class Analyzer(object):
    """docstring for Analyzer"""

    def __init__(self, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        super(Analyzer, self).__init__()
        self.psf = psf
        self.dcd = dcd
        self.cutoff = cutoff
        self.hbondcutoff = hbondcutoff
        self.hbondcutangle = hbondcutangle
        self.sel1text = sel1text
        self.sel2text = sel2text

    def runFrameScan(self):
        self.contactResults = self.analyze_psf_dcd(self.psf, self.dcd, self.cutoff, self.hbondcutoff,
                                                   self.hbondcutangle, self.sel1text, self.sel2text)

    def runContactAnalysis(self,map1, map2):
        self.finalAccumulatedContacts = self.analyze_contactResultsWithMaps(self.contactResults, map1, map2)
        return deepcopy(self.finalAccumulatedContacts)

    def setTrajectoryData(self,resname_array,resid_array,name_array,type_array,segids,backbone):
        self.resname_array = resname_array
        self.resid_array = resid_array
        self.name_array = name_array
        self.type_array = type_array
        self.segids = segids
        self.backbone = backbone

    def getTrajectoryData(self):
        return [self.resname_array,self.resid_array,self.name_array,self.type_array,self.segids,self.backbone]
    # utility for memory measurement
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']

    def humansize(self,nbytes):
        if nbytes == 0: return '0 B'
        i = 0
        while nbytes >= 1024 and i < len(suffixes) - 1:
            nbytes /= 1024.
            i += 1
        f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
        return '%s %s' % (f, suffixes[i])

    # find a string in s between the strings first and last
    def find_between(self,s, first, last):
        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

    # weight function to score contact distances
    def weight_function(self,value):
        return (1.0) / (1.0 + np.exp(5.0 * (value - 4.0)))

    ## maps contain information wether to consider an atom's field for contact accumulation
    # map1 and map2 contain six boolean values each, cf. AccumulationMapIndex
    # for a given contact, the corresponding value to a field is written to keys1 and keys2, respectively
    # example input:
    # map1 = [0,0,0,1,1,0]
    # map2 = [0,0,0,1,1,0], meaning that residue and resname should be used for contact accumulation
    # contact: idx1,idx2
    # results: (example!)
    #	keys1=["none","none","none","14", "VAL", "none"]
    #	keys2=["none","none","none","22", "ILE, "none"]
    def makeKeyArraysFromMaps(self,map1, map2, contact):
        idx1 = contact.idx1
        idx2 = contact.idx2
        counter = 0
        keys1 = []
        for val in map1:
            if val == 1:
                if counter == AccumulationMapIndex.index:
                    keys1.append(idx1)
                elif counter == AccumulationMapIndex.atype:
                    keys1.append(self.type_array[idx1])
                elif counter == AccumulationMapIndex.name:
                    keys1.append(self.name_array[idx1])
                elif counter == AccumulationMapIndex.resid:
                    keys1.append(self.resid_array[idx1])
                elif counter == AccumulationMapIndex.resname:
                    keys1.append(self.resname_array[idx1])
                elif counter == AccumulationMapIndex.segid:
                    keys1.append(self.segids[idx1])
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
                    keys2.append(self.type_array[idx2])
                elif counter == AccumulationMapIndex.name:
                    keys2.append(self.name_array[idx2])
                elif counter == AccumulationMapIndex.resid:
                    keys2.append(self.resid_array[idx2])
                elif counter == AccumulationMapIndex.resname:
                    keys2.append(self.resname_array[idx2])
                elif counter == AccumulationMapIndex.segid:
                    keys2.append(self.segids[idx2])
            else:
                keys2.append("none")
            counter += 1
        return [keys1, keys2]

    # convert a key back to two key arrays
    # cf. comments on makeKeyFromKeyArrays and makeKeyArraysFromMaps
    # "inverse" function of makeKeyFromKeyArrays
    def makeKeyArraysFromKey(self,key):
        keystring1, keystring2 = key.split("-")
        mapping = AccumulationMapIndex.mapping
        maximal = len(mapping)
        key1 = []
        for i in range(0, maximal):
            current = mapping[i]
            if current not in keystring1:
                key1.append("none")
                continue
            if i == (maximal - 1):
                key1.append(keystring1[keystring1.index(current) + len(current):])
                break
            nextCurrent = mapping[i + 1]
            if nextCurrent not in keystring1:
                nxt = ""
                for k in mapping[i + 1:]:
                    if k in keystring1[keystring1.index(current) + len(current):]:
                        nxt = k
                        break
                if nxt != "":
                    key1.append(keystring1[keystring1.index(current) + len(current):keystring1.index(nxt)])
                else:
                    key1.append(keystring1[keystring1.index(current) + len(current):])
                continue
            else:
                currentValue = self.find_between(keystring1, current, nextCurrent)
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
            if i == (maximal - 1):
                key2.append(keystring2[keystring2.index(current) + len(current):])
                break
            nextCurrent = mapping[i + 1]
            if nextCurrent not in keystring2:
                nxt = ""
                for k in mapping[i + 1:]:
                    if k in keystring2[keystring2.index(current) + len(current):]:
                        nxt = k
                        break
                if nxt != "":
                    key2.append(keystring2[keystring2.index(current) + len(current):keystring2.index(nxt)])
                else:
                    key2.append(keystring2[keystring2.index(current) + len(current):])
                continue
            else:
                currentValue = self.find_between(keystring2, current, nextCurrent)
                if currentValue == "":
                    key2.append("none")
                else:
                    key2.append(currentValue)
        return [key1, key2]

    ## input two key arrays as explained above
    #	example:
    #	keys1=["none","none","none","14", "VAL", "none"]
    #	keys2=["none","none","none","22", "ILE, "none"]
    #	returns a human readable key with the mapping identifiers in AccumulationMapIndex
    #	in the given example data:
    #	key="r.14rn.VAL-r.22rn.ILE"
    #	key is used to accumulated contacts in a dictionary (= a contact's unique identifier)
    def makeKeyFromKeyArrays(self,key1, key2):
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

    def analyze_psf_dcd(self,psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        # reading topology/parameter CHARMM files for setting AtomTypes and AtomHBondTypes
        heavyatomlines = []
        heavyatoms = []
        pars = open(os.path.dirname(os.path.abspath(__file__))+'/testpar.prm', 'r')
        for line in pars:
            if re.match("MASS", line):
                heavyatomlines.append(line.rstrip())
        for atomline in heavyatomlines:
            # read new AtomType and its corresponding AtomHBondType from file
            atype = AtomType.parseParameterFileString(atomline)
            heavyatoms.append(atype)

        ### config (GUI settings!)

        # cutoff for contact measurement
        # cutoff = 5.0
        # cutoff for H-A distance, usually 3-3.5?(check in literature again!)
        # hbondcutoff = 3.5
        # cutoff angle for hbond, check literature again
        # hbondcutangle = 120
        ## selection texts (MDAnalysis format, not VMD!)
        # to scan for hydrogen bonds, please do NOT give selections without hydrogen atoms!
        # the computational effort is managed by the algorithm in an appropriate manner!
        # sel1text = "segid RN11"
        # sel2text = "segid UBQ"

        # input psf file (GUI: file picker)
        # psf = "rpn11_ubq_interface-ionized.psf"
        # input dcd file (GUI: file picker)
        # dcd = "short.dcd"
        # dcd = "rpn11_ubq_50ns.dcd"

        ## input maps for contact accumulation
        # boolean values, check AccumulationMapIndex for meaning!
        # map1 = [0,0,0,1,1,0]
        # map2 = [0,0,0,1,1,0]

        ####### MAIN ALGORITHM
        ######################
        # load psf and dcd file in memory
        u = MDAnalysis.Universe(psf, dcd)
        # define selections according to sel1text and sel2text
        sel1 = u.select_atoms(sel1text)
        sel2 = u.select_atoms(sel2text)
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
        self.resname_array = []
        self.resid_array = []
        self.name_array = []
        self.type_array = []
        self.bonds = []
        self.segids = []
        self.backbone = []
        for atom in all_sel.atoms:
            self.resname_array.append(atom.resname)
            self.resid_array.append(atom.resid)
            self.name_array.append(atom.name)
            self.type_array.append(atom.type)
            self.bonds.append(atom.bonds)
            self.segids.append(atom.segid)
        for atom in backbone_sel:
            self.backbone.append(atom.index)
        # show trajectory information and selection information
        print "trajectory with %d frames loaded" % len(u.trajectory)
        print len(sel1.coordinates()), len(sel2.coordinates())
        contactResults = []
        # loop over trajectory
        self.totalFrameNumber = len(u.trajectory)
        for ts in u.trajectory:
            currentFrameContacts = []
            frame = ts.frame
            self.currentFrameNumber = ts.frame
            print frame
            result = np.ndarray(shape=(len(sel1.coordinates()), len(sel2.coordinates())), dtype=float)
            # distarray is the distance matrix between all atoms in sel1 and sel2
            # row = sel1, column = sel2
            distarray = distances.distance_array(sel1.coordinates(), sel2.coordinates(), box=None, result=result)
            contacts = np.where(distarray <= cutoff)
            # idx1 and idx2 correspond to a row,column in contacts, respectively
            # they do NOT correspond to a global atom index!
            for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
                convindex1 = indices1[idx1]  # idx1 converted to global atom indexing
                convindex2 = indices2[idx2]  # idx2 converted to global atom indexing
                # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered, hydrogen bonds can still be detected!
                if re.match("H(.*)", self.name_array[convindex1]) or re.match("H(.*)", self.name_array[convindex2]):
                    continue
                    # distance between atom1 and atom2
                distance = distarray[idx1, idx2]
                weight = self.weight_function(distance)
                # read AtomHBondType from heavyatoms list
                type1 = next((x.htype for x in heavyatoms if x.name == self.type_array[convindex1]), AtomHBondType.none)
                type2 = next((x.htype for x in heavyatoms if x.name == self.type_array[convindex2]), AtomHBondType.none)
                ## HydrogenBondAlgorithm
                # TODO: outsource to another file?
                # elongates the current loop a lot...
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
                        # search for hatom, check numbering in bond!!!!!!!!!!
                        b1 = self.bonds[convindex1]
                        b2 = self.bonds[convindex2]
                        # search for hydrogen atoms bound to atom 1
                        bondcount1 = 0
                        hydrogenAtomsBoundToAtom1 = []
                        for b in b1.types():
                            hydrogen = next((x for x in b if x.startswith("H")), 0)
                            # print b
                            if hydrogen != 0:
                                # print "h bond to atom1"
                                bondindices1 = b1.to_indices()[bondcount1]
                                hydrogenidx = next(
                                    ((j + 1) for j in bondindices1 if self.type_array[j + 1].startswith("H")), -1)
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
                                hydrogenidx = next(
                                    ((k + 1) for k in bondindices2 if self.type_array[k + 1].startswith("H")), -1)
                                if hydrogenidx != -1:
                                    # print type_array[hydrogenidx]
                                    hydrogenAtomsBoundToAtom2.append(hydrogenidx)
                            bondcount2 += 1
                        # check hbond criteria for hydrogen atoms bound to first atom
                        for global_hatom in hydrogenAtomsBoundToAtom1:
                            conv_hatom = indices1.index(global_hatom)
                            typeHeavy = next((x.htype for x in heavyatoms if x.name == self.type_array[convindex2]),
                                             AtomHBondType.none)
                            if typeHeavy == AtomHBondType.acc and (distarray[conv_hatom, idx2] <= hbondcutoff):
                                donorPosition = sel1.coordinates()[idx1]
                                hydrogenPosition = sel1.coordinates()[conv_hatom]
                                acceptorPosition = sel2.coordinates()[idx2]
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
                            typeHeavy = next((x.htype for x in heavyatoms if x.name == self.type_array[convindex1]),
                                             AtomHBondType.none)
                            if typeHeavy == AtomHBondType.acc and (distarray[idx1, conv_hatom] <= hbondcutoff):
                                donorPosition = sel2.coordinates()[idx2]
                                hydrogenPosition = sel2.coordinates()[conv_hatom]
                                acceptorPosition = sel1.coordinates()[idx1]
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
                                    # print str(convindex1) + " " + str(convindex2)
                                    # print "hbond found: %d,%d,%d"%(convindex2,global_hatom,convindex1)
                                    # print angle
                                    # finalize
                newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1),
                                             int(convindex2),
                                             hydrogenBonds)
                currentFrameContacts.append(newAtomContact)
            contactResults.append(currentFrameContacts)
        return contactResults

    def analyze_contactResultsWithMaps(self,contactResults, map1, map2):
        #################################################
        ######## contactResults evaluation
        # only depending on map1, map2
        # part can be run without running the contact analysis algorithm again, as it just prepares the results for displaying
        #################################################

        ### Data structure to process:
        # contactResults (list of frames)
        # ---> frame (list of AtomContacts)
        # --------> AtomContact

        frame_contacts_accumulated = []
        # frame_contacts_accumulated (list of frames)
        # ---> frame_dict (dict)
        # --------> key vs. TempContactAccumulate

        # list of all contacts keys (= unique identifiers, determined by the given maps)
        start = time.time()
        allkeys = []
        for frame in contactResults:
            currentFrameAcc = {}
            for cont in frame:
                key1, key2 = self.makeKeyArraysFromMaps(map1, map2, cont)
                key = self.makeKeyFromKeyArrays(key1, key2)
                if key in currentFrameAcc:
                    currentFrameAcc[key].fscore += cont.weight
                    currentFrameAcc[key].contributingAtomContacts.append(cont)
                    if cont.idx1 in self.backbone:
                        currentFrameAcc[key].bb1score += cont.weight
                    else:
                        currentFrameAcc[key].sc1score += cont.weight
                    if cont.idx2 in self.backbone:
                        currentFrameAcc[key].bb2score += cont.weight
                    else:
                        currentFrameAcc[key].sc2score += cont.weight
                else:
                    currentFrameAcc[key] = TempContactAccumulate(key1, key2)
                    currentFrameAcc[key].fscore += cont.weight
                    currentFrameAcc[key].contributingAtomContacts.append(cont)
                    if cont.idx1 in self.backbone:
                        currentFrameAcc[key].bb1score += cont.weight
                    else:
                        currentFrameAcc[key].sc1score += cont.weight
                    if cont.idx2 in self.backbone:
                        currentFrameAcc[key].bb2score += cont.weight
                    else:
                        currentFrameAcc[key].sc2score += cont.weight
                if not key in allkeys:
                    allkeys.append(key)
            frame_contacts_accumulated.append(currentFrameAcc)
        accumulatedContactsDict = {}
        stop = time.time()
        print stop - start
        # accumulatedContactsDict (dict)
        # ---> key vs. list of TempContactAccumulated
        #
        # loop fills gaps with zero-score TempContactAccumulate of key if key is not occuring in a frame
        # provides clean data!
        start = time.time()
        for key in allkeys:
            accumulatedContactsDict[key] = []
            for frame_dict in frame_contacts_accumulated:
                if not key in frame_dict:  # puts empty score TempContactAccumulate in dict
                    key1, key2 = self.makeKeyArraysFromKey(key)
                    emptyCont = TempContactAccumulate(key1, key2)
                    emptyCont.fscore = 0
                    frame_dict[key] = emptyCont
                accumulatedContactsDict[key].append(frame_dict[key])

                # make a list of AccumulatedContacts from accumulatedContactsDict
        # probably, there is a much easier way to do that, but I am too tired at the moment and it works, though... (M)
        finalAccumulatedContacts = []  # list of AccumulatedContacts
        for key in accumulatedContactsDict:
            key1, key2 = self.makeKeyArraysFromKey(key)
            acc = AccumulatedContact(key1, key2)
            for tempContact in accumulatedContactsDict[key]:
                acc.addScore(tempContact.fscore)
                acc.addContributingAtoms(tempContact.contributingAtomContacts)
                acc.bb1 += tempContact.bb1score
                acc.bb2 += tempContact.bb2score
                acc.sc1 += tempContact.sc1score
                acc.sc2 += tempContact.sc2score
            finalAccumulatedContacts.append(acc)
            print key, acc.bb1, acc.bb2, acc.sc1, acc.sc2
            print len(acc.scoreArray)
        stop = time.time()
        print stop - start
        return finalAccumulatedContacts
