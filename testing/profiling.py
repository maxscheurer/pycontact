from PyContact.core.ContactAnalyzer import *
from PyContact.exampleData.datafiles import DCD, PSF
import MDAnalysis as mda
dcdfile = DCD
psffile = PSF
analyzer = Analyzer(psffile, dcdfile, 5.0, 2.5, 120, "segid RN11", "segid UBQ")
analyzer.runFrameScan()
