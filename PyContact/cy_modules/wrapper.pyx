cdef extern from "src/manager.hh":
    cdef cppclass C_GPUAdder "GPUAdder":
        C_GPUAdder()

cdef class GPUAdder:
    cdef C_GPUAdder* g

    def __cinit__(self):
        self.g = new C_GPUAdder()
