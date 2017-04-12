from cython.view cimport array as cvarray
import numpy as np

cdef extern from "src/gridsearch.C":
  double sasa_grid(const float *pos, int natoms, float pairdist, int allow_double_counting, int maxpairs,
                   const float *radius, const int npts, double srad, int pointstyle,
                   int restricted, const int* restrictedList)

cdef extern from "src/gridsearch.C":
  int* find_within(const float *xyz, int *flgs, int *others, int num, float r)

def cy_sasa(npcoords, natoms, pairdist, allow_double_counting, maxsize, nprad,
  surfacePoints, probeRadius, pointstyle,
  restricted,restrictedList):
  cdef int [::1] cy_restrictedList = restrictedList
  cdef float [::1] cy_radius = nprad
  cdef float [::1] c_coords = npcoords[0]
  cdef float sasa = sasa_grid(&c_coords[0], natoms, pairdist, 0, -1, &cy_radius[0], surfacePoints, probeRadius,
                              pointstyle, restricted, &cy_restrictedList[0])
  return sasa

def cy_find_within(xyz, flgs, others, num, r):
  cdef float [::1] cy_pos = xyz[0]
  cdef int [::1] cy_flags = flgs
  cdef int [::1] cy_others = others
  cdef int[::1] result = <int[:num]> find_within(&cy_pos[0],&cy_flags[0], &cy_others[0], num, r)
  return np.asarray(result)
