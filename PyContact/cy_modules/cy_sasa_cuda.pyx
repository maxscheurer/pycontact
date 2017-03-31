cdef extern from "src/sasaCuda.hh":
    double calculate_sasa_cuda(const double test)


def cy_sasa_cuda(test):
  cdef double sasa = calculate_sasa_cuda(test)
  return sasa
