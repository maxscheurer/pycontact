import multiprocessing

from .ContactAnalyzer import ContactAnalyzer
from .ContactTrajectories import *

class ContactManager:
    """
    Main manager for all contact data and trajectories


    """
    def __init__(self, topology, trajectories, cutoff, hbondcutoff,
                 hbondcutangle, sel1text, sel2text):
        self.topology = topology
        if not isinstance(trajectories, list):
            raise TypeError("Trajectories must be list.")
        self.trajectories = trajectories
        self.trajectoryScanParameters = TrajectoryScanParameters(cutoff,
                                                                 hbondcutoff,
                                                                 hbondcutangle,
                                                                 sel1text,
                                                                 sel2text)

        self.atomicContactTrajectories = []
        self.accumulatedContactTrajectories = []

    def readTrajectories(self, nproc):
        for trajectory in self.trajectories:
            atomicContactTrajectory = ContactAnalyzer.runFrameScan(self.topology,
                                                                   trajectory,
                                                                   self.trajectoryScanParameters,
                                                                   nproc)
            self.atomicContactTrajectories.append(atomicContactTrajectory)

    def accumulateContacts(self, map1, map2):
        for atomicTrajectory in self.atomicContactTrajectories:
            ct = ContactAnalyzer(atomicTrajectory, map1, map2).accumulateContacts()
            self.accumulatedContactTrajectories.append(ct)

class TrajectoryScanParameters:
    def __init__(self, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        self.cutoff = cutoff
        self.hbondcutoff = hbondcutoff
        self.hbondcutangle = hbondcutangle
        self.sel1text = sel1text
        self.sel2text = sel2text
