import collections
import numpy as np
compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

class ResidueType:
    positive, negative, nonpolar, polar, other = range(5)
    mapping = {"asp": negative, "arg": positive, "lys": positive, "glu": negative}

class ContactType:
    saltbr, hydrophobic, other = range(3)
    mapping = [[ResidueType.positive, ResidueType.negative],[ResidueType.nonpolar, ResidueType.nonpolar]]
    colors = ["rgb(255, 0,0)", "rgb(0, 0,255)", "rgb(255, 255 ,255)"]


class BackboneSidechainType:
    contactsBb, contactsSc = range(2)


class BackboneSidechainContactType:
    bb_only, sc_only, both = range(3)
    mapping = [[BackboneSidechainType.contactsBb, BackboneSidechainType.contactsBb],[BackboneSidechainType.contactsBb, BackboneSidechainType.contactsSc], [BackboneSidechainType.contactsSc, BackboneSidechainType.contactsSc]]
    colors = [[0,200,200],[200,200,0],[0,200,0]]


class Residue:
    def __init__(self, name, bb, sc):
        self.name = name.lower()
        self.bb = float(bb)
        self.sc = float(sc)
        if self.bb > self.sc:
            self.contactsBy = BackboneSidechainType.contactsBb
        else:
            self.contactsBy = BackboneSidechainType.contactsSc
        # self.bbScRatio = self.bb/self.sc
        if self.name in ResidueType.mapping:
            self.type = ResidueType.mapping[self.name]
        else:
            self.type = ResidueType.other

class Contact:
    def __init__(self, resA, residA, resB, residB, bb1, sc1, bb2, sc2, scoreArray):
        self.resA = resA
        self.resB = resB
        self.residA = residA
        self.residB = residB
        self.scoreArray = scoreArray
        self.title = self.resA + self.residA + " - " + self.resB + self.residB
        self.residueA = Residue(self.resA,bb1,sc1)
        self.residueB = Residue(self.resB,bb2,sc2)
        self.type = determine_ctype(self.residueA, self.residueB)
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

def determine_ctype(resA, resB):
    if compare([resA.type,resB.type],ContactType.mapping[ContactType.saltbr]):
        return ContactType.saltbr
    else:
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