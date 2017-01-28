'''
    Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
'''
from __future__ import print_function
import collections

import numpy as np
from PyQt5.QtGui import QColor

from ..db.DbReader import dict_factory, read_residue_db


compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

# NOTE: needed for sasa calculation, we should add more...
# CHARMM radii, from VMD
vdwRadii = {"H": 1.0,
            "C": 1.5,
            "N": 1.399999976158142,
            "O": 1.2999999523162842,
            "F": 1.47,
            "Mg": 1.73,
            "P": 1.8,
            "S": 1.899999976158142}

def vdwRadius(atomType):
        return vdwRadii.get(atomType, 1.5)


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


# contains all information of a contact, mapped by key1/key2
class AccumulatedContact(object):
    """Contains information about a contact accumulated from AtomContact to display in GUI"""

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
        """returns the title of the AccumulatedContact to be displayed in contact's label"""
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
            if contactActive is False and score > threshold:
                contactActive = True
                contactTime += ns_per_frame
            elif contactActive is True and score > threshold:
                contactTime += ns_per_frame
            elif contactActive is True and score <= threshold:
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

    def contactTypeAsShortcut(self):
        return ContactType.shortcut[self.determine_ctype()]

    def determine_ctype(self):
        #only works if both maps contain resname
        r1 = self.key1[AccumulationMapIndex.resname].lower()
        r2 = self.key2[AccumulationMapIndex.resname].lower()
        # if r1 == "none" or r2 == "none":
            # return ContactType.other
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
        ishbond = 0
        hb = self.hbondFramesScan()
        for bla in hb:
            if bla > 0:
                ishbond = 1
                break

        if self.atom1contactsBy == BackboneSidechainType.contactsSc and self.atom2contactsBy == BackboneSidechainType.contactsSc:
            # check for saltbridge
            if (scpol1 == SideChainPolarity.positive and scpol2 == SideChainPolarity.negative) or \
                    (scpol2 == SideChainPolarity.positive and scpol1 == SideChainPolarity.negative):
                return ContactType.saltbr
            # check for hydrophobic contact
            if scpol1 == SideChainPolarity.nonpolar and scpol2 == SideChainPolarity.nonpolar:
                if ishbond == 0:
                    return ContactType.hydrophobic
            if ishbond == 1:
                return ContactType.hbond


            return ContactType.other

        if ishbond == 1:
            return ContactType.hbond
        else:
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
        print("frame: %d, dist: %f, weight: %f, idx1: %d, idx2: %d" % (
            self.frame, self.distance, self.weight, self.idx1, self.idx2))


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
        print("donor: %d, acceptor: %d, hydrogen: %d, dist: %f, angle: %f, usedCutoffDist: %f, usedCutoffAngle: %f" % (
            self.donorIndex, self.acceptorIndex, self.hydrogenIndex, self.acceptorHydrogenDistance, self.angle,
            self.usedCutoffDist, self.usedCutoffAngle))

## enum and mapping for atom properties
# used to dynamically define keys and have a bijective nomenclature for all properties
class AccumulationMapIndex():
    index, atype, name, resid, resname, segid = range(6)
    mapping = ["i.", "t.", "nm.", "r.", "rn.", "s."]
    vmdsel = ["index", "type", "name", "resid", "resname", "segname"]

class ContactType:
    saltbr, hydrophobic, hbond, other = range(4)
    shortcut = ["saltbr","hydrophobic","hbond","other"]
    colors = ["rgba(255, 0,0, 50)", "rgba(0, 0,255, 50)", "rgba(255, 0 ,255, 50)", "rgba(255, 255 ,255, 50)"]
    qcolors = [QColor(255, 0, 0, 50), QColor(0, 0, 255, 50), QColor(255, 0, 255, 50), QColor(255, 255, 255, 50)]

class SideChainPolarity:
    nonpolar, positive, negative, polar, other = range(5)
    mapping = {"nonpolar": nonpolar, "positive": positive, "negative": negative, "polar": polar, "other": other}

class BackboneSidechainType:
    contactsBb, contactsSc = range(2)


class BackboneSidechainContactType:
    bb_only, both, sc_only = range(3)
    mapping = [[BackboneSidechainType.contactsBb, BackboneSidechainType.contactsBb],[BackboneSidechainType.contactsBb, BackboneSidechainType.contactsSc], [BackboneSidechainType.contactsSc, BackboneSidechainType.contactsSc]]
    colors = [[0,200,200],[200,200,0],[0,200,0]]

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
