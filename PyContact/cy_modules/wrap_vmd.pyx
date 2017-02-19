cdef extern from "vmd_wrapper.cpp":
  void start_vmd(int pt, char* vmdcommand)
  void stop_vmd()
  void send(char * command)


def start(port, command):
  start_vmd(port, command)

def stop():
  stop_vmd()

def send_command(command):
  send(command)
