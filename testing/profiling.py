#python2 -m cProfile -o pycontact.profile -s cumtime profiling.py
#visualize output with snakeviz
from PyContact.core.ContactAnalyzer import *
from PyContact.core.aroundPatch import AroundSelection
from PyContact.exampleData.datafiles import DCD, PSF
import MDAnalysis as mda
dcdfile = DCD
psffile = PSF
analyzer = Analyzer(psffile, dcdfile, 5.0, 2.5, 120, "segid RN11", "segid UBQ")
analyzer.runFrameScan()
