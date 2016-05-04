import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD   # test trajectory
import numpy.linalg

u = MDAnalysis.Universe(PSF,DCD)  # always start with a Universe
nterm = u.s4AKE.N[0]  # can access via segid (s4AKE) and atom name
cterm = u.s4AKE.C[-1]             # ... takes the last atom named 'C'
bb = u.selectAtoms('protein and backbone')  # a selection (AtomGroup)

for ts in u.trajectory:     # iterate through all frames
    r = cterm.pos - nterm.pos # end-to-end vector from atom positions
    d = numpy.linalg.norm(r)  # end-to-end distance
    rgyr = bb.radiusOfGyration()  # method of AtomGroup
    print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
          ts.frame, d, rgyr))