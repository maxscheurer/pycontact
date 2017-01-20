from cython.view cimport array as cvarray


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


cdef extern from "gridsearch.C":
  int test_function(int i)


cdef extern from "gridsearch.C":
  double sasa_grid(const float *pos,int natoms, float pairdist, int allow_double_counting, int maxpairs, const float *radius,const int npts, double srad, int pointstyle, int restricted,  const int* restrictedList)

def bla(int i):
  return test_function(i)


def test_sasa():
  import MDAnalysis
  import numpy as np
  psf = "/mnt/workspace/pycontactData/nowater.psf"
  dcd = "/mnt/workspace/pycontactData/trajectory_short.dcd"

  # load psf and trajectory, make lists with radii and coordinates
  u = MDAnalysis.Universe(psf, dcd)

  probeRadius = 1.4

  # seltext = "segid UBQ"
  # resseltext = "segid UBQ and same residue as around 5.0 (segid RN11)"
  seltext = "segid COH3"
  seltext2 = "segid DOC3"
  resseltext = "segid COH3 and around 5 segid DOC3"
  perres = 0

  # 0=spiral, 1=random (VMD)
  pointstyle = 1
  # number of points to approximate the sphere
  surfacePoints = 50
  # pair distance
  pairdist = 2 * (2.0 + 1.4)

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
  nprad = np.array(radius, dtype=np.float32)
  restrictedList = np.array(restrictedList, dtype=np.int32)
  cdef int [::1] cy_restrictedList = restrictedList
  cdef float [::1] cy_radius = nprad

  input_coords = []
  ressel = u.select_atoms(resseltext)
  for ts in u.trajectory:
      # print("restricted: ", len(ressel.atoms))
      input_coords.append(selection.positions)

  cdef float [::1] c_coords;
  result = []
  for c in input_coords:
        coords = np.reshape(c, (1, natoms * 3))
        npcoords = np.array(coords, dtype=np.float32)
        print(npcoords[0])
        c_coords = npcoords[0]
        sasa = sasa_grid(&c_coords[0], natoms, pairdist, 0, -1, &cy_radius[0] ,surfacePoints, probeRadius, pointstyle, restricted, &cy_restrictedList[0])
        result.append(sasa)
  return result
