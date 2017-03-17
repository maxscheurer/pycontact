"""contains classes in order to faciliate scripting the PyContact package."""
from PyContact.core.ContactAnalyzer import *
from PyContact.core.aroundPatch import AroundSelection
from PyContact.core.DataHandler import DataHandler


class JobConfig:
    """Configuration/Settings for a PyContact contact analysis job"""
    def __init__(self, cutoff, hbondcutoff, hbondcutangle, map1, map2, sel1, sel2):
        self.cutoff = cutoff
        self.hbondcutoff = hbondcutoff
        self.hbondcutangle = hbondcutangle
        self.map1 = map1
        self.map2 = map2
        self.sel1 = sel1
        self.sel2 = sel2


class PyContactJob:
    """Job class that is given a JobConfig and input files and handles
    running/analyzing the job as such."""
    def __init__(self, topo, traj, name, configuration):
        self.topo = topo
        self.traj = traj
        self.name = name
        self.configuration = configuration
        self.analyzer = None

    def runJob(self, ncores=1):
        """Runs contact analysis and accumulation with ncores threads."""
        print("Running job: " + self.name + " on " + str(ncores) + " cores.")
        self.analyzer = Analyzer(self.topo,self.traj, self.configuration.cutoff, self.configuration.hbondcutoff, self.configuration.hbondcutangle, self.configuration.sel1, self.configuration.sel2)
        self.analyzer.runFrameScan(ncores)
        self.analyzer.runContactAnalysis(self.configuration.map1, self.configuration.map2, ncores)

    def writeSessionToFile(self, fname=""):
        """Writes the current analysis session to a file with either self.name + .session or fname as output filename"""
        if fname != "":
            DataHandler.writeSessionToFile(fname, self.analyzer)
        else:
            DataHandler.writeSessionToFile(self.name + ".session", self.analyzer)
        print("Wrote session to file")
