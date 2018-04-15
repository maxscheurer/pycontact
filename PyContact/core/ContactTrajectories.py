import numpy as np

from .Biochemistry import (BackboneSidechainType,
                           BackboneSidechainContactType,
                           makeHumanReadableTitle)
from .KeyManager import KeyManager

class AtomicContactTrajectory:
    """docstring for [object Object]."""
    def __init__(self, contacts, resname_array, resid_array, name_array,
                 segids, backbone):
        self.contacts = contacts
        self.resname_array = resname_array
        self.resid_array = resid_array
        self.name_array = name_array
        self.segids = segids
        self.backbone = backbone



class AccumulatedContactTrajectory:
    def __init__(self, keys, contactScores, bbScores, scScores,
                 hbonds):
        """
        Data container for trajectory contact data

        Parameters
        ----------
        keys: np.ndarray
            array of contact keys
        contactScores: np.ndarray
            accumulated contact scores per frame, shape(contact,frame)
        bbScores: tuple of 2 np.ndarrays
            backbone scores
        scScores: tuple of 2 np.ndarrays
            side chain scores
        hbonds: np.ndarray
            hydrogen bonds, shape(contact,frame)

        """
        self.keys = keys
        self.contactScores = contactScores
        self.bbScores1, self.bbScores2 = bbScores
        self.scScores1, self.scScores2 = scScores
        self.hbonds = hbonds
        self.numberOfFrames = len(self.keys)
        self.backboneSideChainTypes = np.zeros(len(self.contactScores), dtype=np.int8)
        self.titles = np.array([])
        self.determineBackboneSidechainType()
        self.makeTitles()

    def determineBackboneSidechainType(self):
        """Sets the Backbone-Sidechain type."""
        atom1 = self.bbScores1 > self.scScores1
        atom2 = self.bbScores2 > self.scScores2
        bb = (atom1 == atom2)

        self.backboneSideChainTypes[np.where(bb == False)] = BackboneSidechainContactType.both
        self.backboneSideChainTypes[np.where((atom1 == True) & (atom2 == True))] = BackboneSidechainContactType.bb_only
        self.backboneSideChainTypes[np.where((atom1 == False) & (atom2 == False))] = BackboneSidechainContactType.sc_only

    def makeTitles(self):
        self.titles = np.array([])
        for k in self.keys:
            key_array = KeyManager.makeKeyArraysFromKey(k)
            self.titles = np.append(self.titles,
                                    makeHumanReadableTitle(*key_array))
