import MDAnalysis as md
from MDAnalysis.lib.formats.libdcd import DCDFile

psf = "../rpn11_ubq_interface-ionized.psf"
dcd = "../short.dcd"

# u = md.Universe(psf,dcd)
dcdfile = DCDFile(dcd)
