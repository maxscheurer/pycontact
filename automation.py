import pickle

from PyContact.core.ContactAnalyzer import *
from PyContact.core.aroundPatch import AroundSelection
from PyContact.core.multi_trajectory import run_load_parallel

job = ["/run/media/max/max_data/dodecin/MD_all_woMG/run1/dodecin_womg-ionized.psf","/run/media/max/max_data/dodecin/MD_all_woMG/run1/dodecin_womg-unw-al.dcd","ligand_water"]

psffile = job[0]
dcdfile = job[1]
name = job[2]

analyzer = Analyzer(psffile, dcdfile, 5.0, 2.5, 120, "segid LIG", "resname TIP3")
analyzer.runFrameScan()
map1 = [0, 0, 0, 1, 0, 0]
map2 = [0, 0, 0, 0, 1, 0]
contacts = analyzer.runContactAnalysis(map1, map2)

analyzerArgs = [analyzer.psf, analyzer.dcd, analyzer.cutoff, analyzer.hbondcutoff, analyzer.hbondcutangle, analyzer.sel1text, analyzer.sel2text, analyzer.contactResults]
trajArgs = analyzer.getTrajectoryData()
exportDict = {"contacts":contacts,"analyzer":analyzerArgs,"trajectory":trajArgs,"maps":[analyzer.lastMap1,analyzer.lastMap2]}
pickle.dump(exportDict, open(name + ".session", "wb"))
