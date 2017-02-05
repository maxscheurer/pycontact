#include <cmath>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <sys/stat.h>
#include <cstring>
#include <vmdstream/vmdstream.h>

using namespace std;
vmdsock_t vmdsock;


static void start_vmd(int pt) {
  int port = pt;
  vmdsock = newvmdsock("vmd", port);
  vmdstream vmdscript(vmdsock);

  //draw a sphere
  vmdscript << "draw sphere { 0 0 0 } radius 0.5" << endl;

  //render to a file
  // vmdscript << "render snapshot sphere.tga" << endl;

  //quit vmd and exit
  //
  //
  //
}

static void stop_vmd() {
  // vmdscript << "quit" << endl;
  // vmdscript.flush();
  closevmdsock(vmdsock);
}
