import numpy as np

cdef extern from "src/sasaCuda.hh":
    double calculate_sasa_cuda(float *pos, int natoms, float pairdist, float *radius, int npts,
                               float srad, int pointstyle, int restricted, int* restrictedList)


def cy_sasa_cuda(npcoords, natoms, pairdist, nprad, surfacePoints, probeRadius, pointstyle, restricted, restrictedList):

    cdef int [::1] cy_restrictedList = restrictedList
    cdef float [::1] cy_radius = nprad
    cdef float [::1] c_coords = npcoords
    cdef float sasa = calculate_sasa_cuda(&c_coords[0], natoms, pairdist, &cy_radius[0], surfacePoints,
                                          probeRadius, pointstyle, restricted, &cy_restrictedList[0])
    return sasa
