import pickle

from PyContact.core.ContactAnalyzer import *
from PyContact.core.aroundPatch import AroundSelection
from PyContact.core.multi_trajectory import run_load_parallel

class JobConfig:
    def __init__(self, cutoff, hbondcutoff, hbondcutangle, map1, map2, sel1, sel2):
        self.cutoff = cutoff
        self.hbondcutoff = hbondcutoff
        self.hbondcutangle = hbondcutangle
        self.map1 = map1
        self.map2 = map2
        self.sel1 = sel1
        self.sel2 = sel2

class PyContactJob:
    def __init__(self, topo, traj, name, configuration):
        self.topo = topo
        self.traj = traj
        self.name = name
        self.configuration = configuration

    def runJob(self,ncores=1):
        analyzer = Analyzer(self.psffile,self.dcdfile, self.configuration.cutoff, self.configuration.hbondcutoff, self.configuration.hbondcutangle, self.configuration.sel1, self.configuration.sel2)
        analyzer.runFrameScan(ncores)
        analyzer.runContactAnalysis(self.configuration.map1, self.configuration.map2)

    def writeSessionToFile(self):
        # write sth. in DataHandler, pass it analyzer and then write session to file
        pass

jobs = [["/run/media/max/max_data/dodecin/MD_all_woMG/run1/dodecin_womg-ionized.psf","/run/media/max/max_data/dodecin/MD_all_woMG/run1/dodecin_womg-unw-al.dcd","ligand_water"]


analyzerArgs = [analyzer.psf, analyzer.dcd, analyzer.cutoff, analyzer.hbondcutoff, analyzer.hbondcutangle, analyzer.sel1text, analyzer.sel2text, analyzer.contactResults]
trajArgs = analyzer.getTrajectoryData()
exportDict = {"contacts":contacts,"analyzer":analyzerArgs,"trajectory":trajArgs,"maps":[analyzer.lastMap1,analyzer.lastMap2]}
pickle.dump(exportDict, open(name + ".session", "wb"))
