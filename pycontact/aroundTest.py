from __future__ import print_function
import MDAnalysis
from aroundPatch import *
MDAnalysis.core.flags['use_periodic_selections'] = False
MDAnalysis.core.flags['use_KDTree_routines'] = True
psf = "rpn11_ubq_interface-ionized.psf"
pdb = "rpn11_ubq_interface-ionized.pdb"
dcd = "/home/max/Projects/pycontact/short.dcd"

# load psf and trajectory, make lists with radii and coordinates
u = MDAnalysis.Universe(psf, dcd)

seltext="segid UBQ and around 10 segid RN11 or (resid 70-100)"
# selection = u.select_atoms(seltext)
import time
start = time.time()
for ts in u.trajectory:
    selection = u.select_atoms(seltext)
    print("atoms in selection: ", selection.positions.size)
stop = time.time()
print(stop - start)
