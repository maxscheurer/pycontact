#include <cmath>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <sys/stat.h>
#include <cstring>
#include "vmdstream.h"

using namespace std;
vmdsock_t vmdsock;


static void start_vmd(int pt, char* vmdcommand) {
  int port = pt;
  vmdsock = newvmdsock(vmdcommand, port);
  vmdstream vmdscript(vmdsock);
}

static void stop_vmd() {
  vmdstream vmdscript(vmdsock);
  vmdscript << "quit" << endl;
  vmdscript.flush();
  closevmdsock(vmdsock);
}

static void send(char * command) {
  cout << command << endl;
  vmdstream vmdscript(vmdsock);
  vmdscript << command << endl;
  vmdscript.flush();
  vmdscript << "pwd" << endl;
  vmdscript.flush();
}
