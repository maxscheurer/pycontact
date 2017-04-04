import numpy as np

cdef extern from "src/sasaCuda.hh":
    double calculate_sasa_cuda(const float *pos, int natoms, float pairdist, const float *radius, const int npts,
                               double srad, int pointstyle, int restricted, const int* restrictedList)


def cy_sasa_cuda(npcoords, natoms, pairdist, nprad, surfacePoints, probeRadius, pointstyle, restricted, restrictedList):

    cdef int [::1] cy_restrictedList = restrictedList
    cdef float [::1] cy_radius = nprad
    cdef float [::1] c_coords = npcoords[0]
    cdef float sasa = calculate_sasa_cuda(&c_coords[0], natoms, pairdist, &cy_radius[0], surfacePoints,
                                          probeRadius, pointstyle, restricted, &cy_restrictedList[0])
    return sasa
