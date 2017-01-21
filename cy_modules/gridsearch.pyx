from cython.view cimport array as cvarray

cdef extern from "src/gridsearch.C":
  double sasa_grid(const float *pos,int natoms, float pairdist, int allow_double_counting, int maxpairs, const float *radius,const int npts, double srad, int pointstyle, int restricted,  const int* restrictedList)

def cy_sasa(npcoords, natoms, pairdist, allow_double_counting, maxsize, nprad,
  surfacePoints, probeRadius, pointstyle,
  restricted,restrictedList):
  print(probeRadius)
  cdef int [::1] cy_restrictedList = restrictedList
  cdef float [::1] cy_radius = nprad
  cdef float [::1] c_coords = npcoords[0]
  cdef float sasa = sasa_grid(&c_coords[0], natoms, pairdist, 0, -1, &cy_radius[0] ,surfacePoints, probeRadius, pointstyle, restricted, &cy_restrictedList[0])
  return sasa
