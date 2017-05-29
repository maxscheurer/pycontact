cdef extern from "vmd_wrapper.cpp":
  void start_vmd(int pt)
  void stop_vmd()



def start():
  start_vmd(5500)

def stop():
  stop_vmd()
