#usage: mpirun -np 4 python mpi_first.py
from mpi4py import MPI
import numpy as np
import MDAnalysis
from mdanalysis import *
import math
import time
glob_start = time.time()

def chunks(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()



if rank == 0:
	u = MDAnalysis.Universe("rpn11_ubq_interface-ionized.psf", "short.dcd")
	traj_chunks = chunks(u.trajectory,size)
	sel1 = u.select_atoms("segid RN11")
	sel2 = u.select_atoms("segid UBQ")
	indices1 = []
	for at in sel1.atoms:
	    indices1.append(at.index)
	indices2 = []
	for at in sel2.atoms:
	    indices2.append(at.index)
	    # write properties of all atoms to lists
	all_sel = u.select_atoms("all")
	backbone_sel = u.select_atoms("backbone")
	resname_array = []
	resid_array = []
	name_array = []
	type_array = []
	bonds = []
	segids = []
	backbone = []
	for atom in all_sel.atoms:
	    resname_array.append(atom.resname)
	    resid_array.append(atom.resid)
	    name_array.append(atom.name)
	    type_array.append(atom.type)
	    bonds.append(atom.bonds)
	    segids.append(atom.segid)
	for atom in backbone_sel:
	    backbone.append(atom.index)
	comm.send(resname_array,dest=1,tag=12)
else:
	resname_array = comm.recv(tag=12,source=0)


traj = comm.scatter(traj_chunks,root=0)

for ts in traj:
	print ts.frame