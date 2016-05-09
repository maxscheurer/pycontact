'''
    Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
'''
import collections
import numpy as np
from read_db import *
from PyQt5.QtGui import QColor
compare = lambda x, y: collections.Counter(x) == collections.Counter(y)
from analysis_utils import *

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
        self.title = makeKeyFromKeyArrays(key1, key2)
        self.bb1 = 0
        self.sc1 = 0
        self.bb2 = 0
        self.sc2 = 0

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

        if  self.atom1contactsBy == BackboneSidechainType.contactsBb and self.atom2contactsBy == BackboneSidechainType.contactsBb:
            self.bbSidechainType = BackboneSidechainContactType.bb_only
        elif self.atom1contactsBy == BackboneSidechainType.contactsSc and self.atom2contactsBy == BackboneSidechainType.contactsSc:
            self.bbSidechainType = BackboneSidechainContactType.sc_only
        else:
            self.bbSidechainType = BackboneSidechainContactType.both
        return self.bbSidechainType

    def addContributingAtoms(self, contAtoms):
        # append a list of contributing atom to the contributingAtoms list, e.g. when a new frame is added
        self.contributingAtoms.append(contAtoms)  ## used for temporary accumulation of contacts in data analysis

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

class ContactType:
    saltbr, hydrophobic, hbond, other = range(4)
    colors = ["rgba(255, 0,0, 50)", "rgba(0, 0,255, 50)", "rgba(255, 0 ,255, 50)", "rgba(255, 255 ,255, 50)"]
    qcolors = [QColor(255, 0, 0, 50), QColor(0, 0, 255, 50), QColor(255, 0, 255, 50), QColor(255, 255, 255, 50)]

class HBondType:
    none, donor, acceptor, both, other = range(5)
    mapping = {"none": none, "don": donor, "acc": acceptor, "both": both, "other": other}

class SideChainPolarity:
    nonpolar, positive, negative, polar, other = range(5)
    mapping = {"nonpolar": nonpolar, "positive": positive, "negative": negative, "polar": polar, "other": other}

class BackboneSidechainType:
    contactsBb, contactsSc = range(2)


class BackboneSidechainContactType:
    bb_only, both, sc_only = range(3)
    mapping = [[BackboneSidechainType.contactsBb, BackboneSidechainType.contactsBb],[BackboneSidechainType.contactsBb, BackboneSidechainType.contactsSc], [BackboneSidechainType.contactsSc, BackboneSidechainType.contactsSc]]
    colors = [[0,200,200],[200,200,0],[0,200,0]]


class Residue:
    def __init__(self, name, bb, sc, ident):
        self.name = name.lower()
        self.bb = float(bb)
        self.sc = float(sc)
        # identifier = resid in this case
        self.ident = int(ident)
        if self.bb > self.sc:
            self.contactsBy = BackboneSidechainType.contactsBb
        else:
            self.contactsBy = BackboneSidechainType.contactsSc
        # TODO: propably init db first and just call select (speed up)
        #sidechain polarity
        try:
            scpol = str(read_residue_db("scpolarity","name",self.name)[0]["scpolarity"])
            self.scpolarity = SideChainPolarity.mapping[scpol]
        except:
            self.scpolarity = SideChainPolarity.other

        #hydrogen bonds: donor, acceptor, both

        try:
            hbond = str(read_residue_db("hbondtype", "name", self.name)[0]["hbondtype"])
            self.hbondtype = HBondType.mapping[hbond]
        except:
            self.hbondtype = HBondType.other
        # self.printself()

    def printself(self):
        print(self.name + " " + str(self.contactsBy) + " " + str(self.scpolarity) + " " + str(self.hbondtype))

class Contact:
    def __init__(self, resA, residA, resB, residB, bb1, sc1, bb2, sc2, scoreArray):
        self.resA = resA
        self.resB = resB
        self.residA = residA
        self.residB = residB
        self.scoreArray = scoreArray
        self.title = self.resA + self.residA + " - " + self.resB + self.residB
        self.residueA = Residue(self.resA,bb1,sc1,residA)
        self.residueB = Residue(self.resB,bb2,sc2,residB)
        self.contactType = self.determine_ctype()
        self.determineBackboneSidechainType()
        self.mean_score()
        self.median_score()

    def determineBackboneSidechainType(self):
        if compare([self.residueA.contactsBy, self.residueB.contactsBy], BackboneSidechainContactType.mapping[BackboneSidechainContactType.bb_only]):
            self.backboneSideChainType =  BackboneSidechainContactType.bb_only
        elif compare([self.residueA.contactsBy, self.residueB.contactsBy], BackboneSidechainContactType.mapping[BackboneSidechainContactType.sc_only]):
            self.backboneSideChainType = BackboneSidechainContactType.sc_only
        else:
            self.backboneSideChainType = BackboneSidechainContactType.both

    def framenumber(self):
        return len(self.scoreArray)

    def total_time(self, ns_per_frame, threshold):
        time = 0
        for score in self.scoreArray:
            if score > threshold:
                time += ns_per_frame
        self.ttime = time
        return self.ttime



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
            if  i == (len(self.scoreArray) - 1):
                lifeTimes.append(contactTime)
            i += 1
        return lifeTimes

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

    def determine_ctype(self):
        a = self.residueA
        b = self.residueB
        # check if both residues contact by backbone
        if a.contactsBy == BackboneSidechainType.contactsBb and b.contactsBy == BackboneSidechainType.contactsBb:
            return ContactType.hbond

        # check if contact by sidechain and backbone
        if (a.contactsBy == BackboneSidechainType.contactsBb and b.contactsBy == BackboneSidechainType.contactsSc):
            # sc contact by residue B
            if b.hbondtype != HBondType.none:
                return ContactType.hbond
            else:
                return ContactType.other
        elif (a.contactsBy == BackboneSidechainType.contactsSc and b.contactsBy == BackboneSidechainType.contactsBb):
            # sc contact by residue A
            if a.hbondtype != HBondType.none:
                return ContactType.hbond
            else:
                return ContactType.other

        if a.contactsBy == BackboneSidechainType.contactsSc and b.contactsBy == BackboneSidechainType.contactsSc:
            # check for saltbridge
            if (a.scpolarity == SideChainPolarity.positive and b.scpolarity == SideChainPolarity.negative) or \
                (b.scpolarity == SideChainPolarity.positive and a.scpolarity == SideChainPolarity.negative):
                return ContactType.saltbr

            # check for hydrophobic contact
            if a.scpolarity == SideChainPolarity.nonpolar and b.scpolarity == SideChainPolarity.nonpolar:
                if a.hbondtype == HBondType.none or b.hbondtype == HBondType.none:
                    return ContactType.hydrophobic

            # final hbond scan
            if (a.hbondtype == HBondType.donor and b.hbondtype == HBondType.acceptor) or \
                (a.hbondtype == HBondType.acceptor and b.hbondtype == HBondType.donor) or \
                (a.hbondtype == HBondType.both and b.hbondtype == HBondType.both) or \
                (a.hbondtype == HBondType.donor and b.hbondtype == HBondType.both) or \
                (a.hbondtype == HBondType.both and b.hbondtype == HBondType.donor) or \
                (a.hbondtype == HBondType.acceptor and b.hbondtype == HBondType.both) or \
                (a.hbondtype == HBondType.both and b.hbondtype == HBondType.acceptor):
                    return ContactType.hbond

            return ContactType.other

def mean_score_of_contactArray(contacts):
    meanList = []
    for c in contacts:
        meanList = np.concatenate((meanList, c.scoreArray), axis=0)
    return np.mean(meanList)

def median_score_of_contactArray(contacts):
    medianList = []
    for c in contacts:
        medianList = np.concatenate((medianList,c.scoreArray),axis=0)
    return np.median(medianList)