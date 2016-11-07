import MDAnalysis
import numpy as np
from numpy import linalg as la
np.set_printoptions(threshold=np.inf)
from ctypes import cdll
import ctypes
from numpy.ctypeslib import ndpointer
import time
import itertools

# compile: g++ -Wall -g -std=c++11 -fPIC -O -shared gridsearch.C -o libgridsearch.so -lpython2.7

# vdwRadii = {"H": 1.1,
#             "C": 1.7,
#             "N": 1.55,
#             "O": 1.52,
#             "F": 1.47,
#             "Mg": 1.73,
#             "P": 1.8,
#             "S": 1.8}

# CHARMM radii, from VMD
vdwRadii = {"H": 1.0,
            "C": 1.5,
            "N": 1.399999976158142,
            "O": 1.2999999523162842,
            "F": 1.47,
            "Mg": 1.73,
            "P": 1.8,
            "S": 1.899999976158142}

def vdwRadius(atomType):
        return vdwRadii.get(atomType, 1.5)

class SurfaceAnalyser:
    """docstring for SurfaceAnalyser"""

    def computeAsa(self):
        #load shared libraries
        cgrid = cdll.LoadLibrary('./shared/libgridsearch.so')
        search = cgrid.sasa_grid
        search.restype = ctypes.c_double
        search.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"), ctypes.c_int,ctypes.c_float, ctypes.c_int,ctypes.c_int, \
        ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int,ndpointer(ctypes.c_uint32, flags="C_CONTIGUOUS")]

        psf = "rpn11_ubq_interface-ionized.psf"
        pdb = "rpn11_ubq_interface-ionized.pdb"
        dcd = "/home/max/Projects/pycontact/short.dcd"

        # load psf and trajectory, make lists with radii and coordinates
        u = MDAnalysis.Universe(psf, dcd)

        probeRadius = 1.4
        seltext="protein"
        resseltext="segid UBQ and same residue as around 5.0 (segid RN11)"
        perres = 0

        # 0=spiral, 1=random (VMD)
        pointstyle = 1
        # number of points to approximate the sphere
        surfacePoints = 100
        # pair distance
        pairdist = 2*(2.0+1.4)

        if resseltext != "":
            restricted = 1
        else:
            restricted = 0

        selection = u.select_atoms(seltext)

        if perres:
            resids = sorted(set(selection.resids))
            segs = sorted(set(selection.segids))
        else:
            pass

        natoms = len(selection.atoms)
        radius = []
        restrictedList = []
        if restricted:
            ressel = u.select_atoms(resseltext)
            for s in selection.atoms:
                if s in ressel.atoms:
                    restrictedList.append(1)
                else:
                    restrictedList.append(0)
                radius.append(vdwRadius(s.name[0]))
        else:
            for s in selection.atoms:
                radius.append(vdwRadius(s.name[0]))
        natoms = len(selection)
        nprad = np.array(radius,dtype=np.float32)
        restrictedList = np.array(restrictedList,dtype=np.uint32)

        looptime = time.time()
        for frame in u.trajectory:
            coords = np.reshape(selection.positions,(1,natoms*3))
            npcoords = np.array(coords,dtype=np.float32)
            print "start C"
            startC = time.time()
            # sasa_grid(const float *pos,int natoms, float pairdist, int allow_double_counting, int maxpairs, const float *radius,const int npts, double srad, int pointstyle)
            # point style: 0=spiral, 1=random
            asa = search(npcoords, natoms, pairdist, 0,-1,nprad,surfacePoints,probeRadius,pointstyle,restricted,restrictedList)
            stopC = time.time()
            print "time for grid search: ", (stopC-startC)
            print "asa:", asa
        looptimeEnd = time.time()
        print "Loop time: ", looptimeEnd-looptime

a = SurfaceAnalyser()
a.computeAsa()
