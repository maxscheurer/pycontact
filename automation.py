import pickle

from PyContact.core.ContactAnalyzer import *
from PyContact.core.aroundPatch import AroundSelection
from PyContact.core.DataHandler import DataHandler


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
        self.analyzer = None

    def runJob(self,ncores=1):
        print("Running job: " + self.name + " on " + str(ncores) + " cores.")
        self.analyzer = Analyzer(self.topo,self.traj, self.configuration.cutoff, self.configuration.hbondcutoff, self.configuration.hbondcutangle, self.configuration.sel1, self.configuration.sel2)
        self.analyzer.runFrameScan(ncores)
        self.analyzer.runContactAnalysis(self.configuration.map1, self.configuration.map2, ncores)

    def writeSessionToFile(self, fname=""):
        if fname != "":
            DataHandler.writeSessionToFile(fname, self.analyzer)
        else:
            DataHandler.writeSessionToFile(self.name + ".session", self.analyzer)
        print("Wrote session to file")


# jobAllMD = PyContactJob("/run/media/max/max_data/dodecin/MD_all/run1/dodecin_mg-ionized.psf",
#                    "/run/media/max/max_data/dodecin/MD_all/md_all.dcd", "md_all", JobConfig(5.0, 2.5, 120, [0,0,0,0,0,1], [0,0,0,0,0,1], "resname LIG", "resname TIP3"))
#
# jobAllMD.runJob(8)
# jobAllMD.writeSessionToFile()

# jobAllMD_womg = PyContactJob("/run/media/max/max_data/dodecin/MD_all_woMG/run1/dodecin_womg-ionized.psf",
#                    "/run/media/max/max_data/dodecin/MD_all_woMG/md_all_womg.dcd", "md_all_womg", JobConfig(5.0, 2.5, 120, [0,0,0,0,0,1], [0,0,0,0,0,1], "resname LIG", "resname TIP3"))
#
# jobAllMD_womg.runJob(8)
# jobAllMD_womg.writeSessionToFile()

# job_single = PyContactJob("/run/media/max/max_data/dodecin/MD_single_lmf/run1/dodecin_mg_single_lig-ionized.psf",
#                    "/run/media/max/max_data/dodecin/MD_single_lmf/md_single_lmf.dcd", "md_single", JobConfig(5.0, 2.5, 120, [0,0,0,0,0,1], [0,0,0,0,0,1], "resname LMF", "resname TIP3"))
#
# job_single.runJob(8)
# job_single.writeSessionToFile()
#

# jobAllMD_gln = PyContactJob("/run/media/max/max_data/dodecin/MD_all/run1/dodecin_mg-ionized.psf",
#                    "/run/media/max/max_data/dodecin/MD_all/md_all.dcd", "md_all_gln", JobConfig(5.0, 2.5, 120, [0,0,0,0,0,1], [0,0,0,0,0,1], "resname LIG", "resid 55 and protein"))
#
# jobAllMD_gln.runJob(8)
# jobAllMD_gln.writeSessionToFile()


# jobAllMD_womg_gln = PyContactJob("/run/media/max/max_data/dodecin/MD_all_woMG/run1/dodecin_womg-ionized.psf",
#                    "/run/media/max/max_data/dodecin/MD_all_woMG/md_all_womg.dcd", "md_all_womg_gln", JobConfig(5.0, 2.5, 120, [0,0,1,0,0,1], [0,0,0,0,1,0], "resname LIG", "resid 55 and protein"))
#
# jobAllMD_womg_gln.runJob(8)
# jobAllMD_womg_gln.writeSessionToFile()
#
#
# job_single_gln = PyContactJob("/run/media/max/max_data/dodecin/MD_single_lmf/run1/dodecin_mg_single_lig-ionized.psf",
#                    "/run/media/max/max_data/dodecin/MD_single_lmf/md_single_lmf.dcd", "md_single_gln", JobConfig(5.0, 2.5, 120, [0,0,1,0,0,1], [0,0,0,0,1,0], "resname LMF", "resid 55 and protein"))
#
# job_single_gln.runJob(8)
# job_single_gln.writeSessionToFile()
