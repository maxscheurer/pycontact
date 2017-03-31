import numpy as np

cdef extern from "src/sasaCudaKernel.cu":
  double sasa_cuda(const double test)


def cy_sasa(test):
  cdef float sasa = sasa_cuda(test)
  return sasa
