import numpy as np

from .Biochemistry import (BackboneSidechainType,
                           BackboneSidechainContactType,
                           AccumulationMapIndex,
                           ContactType, SideChainPolarity,
                           AminoAcids)
from .KeyManager import KeyManager


def determineContactTypes(keys, hbonds, atom1by, atom2by):
    contactTypes = np.zeros_like(keys, dtype=np.int8)
    idx = 0
    for k, hb in zip(keys, hbonds):
        key1, key2 = KeyManager.makeKeyArraysFromKey(k)
        resname1 = key1[AccumulationMapIndex.resname].lower()
        resname2 = key2[AccumulationMapIndex.resname].lower()

        scpol1 = AminoAcids.scProperties[resname1]
        scpol2 = AminoAcids.scProperties[resname2]

        ishbond = np.any(hb > 0.0)

        contactTypes[idx] = ContactType.other

        if ishbond:
            contactTypes[idx] = ContactType.hbond

        if atom1by[idx] == BackboneSidechainType.contactsSc and atom2by[idx] == BackboneSidechainType.contactsSc:
            if (scpol1 == SideChainPolarity.positive and scpol2 == SideChainPolarity.negative) or \
                    (scpol2 == SideChainPolarity.positive and scpol1 == SideChainPolarity.negative):
                contactTypes[idx] = ContactType.saltbr
            elif scpol1 == SideChainPolarity.nonpolar and scpol2 == SideChainPolarity.nonpolar and ishbond == 0:
                contactTypes[idx] = ContactType.hydrophobic
        idx += 1
    return contactTypes


def makeHumanReadableTitle(key1, key2):
    """returns the title of the AccumulatedContact to be displayed in contact's label"""
    total = []
    for key in [key1, key2]:
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
        tempList = [residueString, atomIndexString, atomNameString, segnameString]
        finishedList = []
        for string in tempList:
            if string != "":
                finishedList.append(string)
        finishedString = " , ".join(finishedList)
        total.append(finishedString)
    return " - ".join(total)


def determineBackboneSidechainTypes(bbScores1, bbScores2,
                                    scScores1, scScores2):
    """Determines the Backbone-Sidechain type."""
    atom1by = bbScores1 > scScores1
    atom2by = bbScores2 > scScores2
    bb = (atom1by == atom2by)

    backboneSideChainTypes = np.zeros_like(bbScores1, dtype=np.int8)

    backboneSideChainTypes[np.where(bb == False)] = BackboneSidechainContactType.both # 1
    backboneSideChainTypes[np.where((atom1by == True) & (atom2by == True))] = BackboneSidechainContactType.bb_only # 0
    backboneSideChainTypes[np.where((atom1by == False) & (atom2by == False))] = BackboneSidechainContactType.sc_only #2
    return (atom1by, atom2by, backboneSideChainTypes)


def makeTitles(keys):
    titles = np.array([])
    for k in keys:
        key_array = KeyManager.makeKeyArraysFromKey(k)
        titles = np.append(titles, makeHumanReadableTitle(*key_array))
    return titles
