import MDAnalysis as md

psf = "../rpn11_ubq_interface-ionized.psf"
dcd = "../short.dcd"

u = md.Universe(psf,dcd)
