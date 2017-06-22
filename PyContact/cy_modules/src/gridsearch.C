#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ResizeArray.h"
#include <vector>

#define FALSE   0
#define TRUE    1
void * find_within_routine( void *v );

struct AtomEntry {
  float x, y, z;
  int index;
  AtomEntry() {}
  AtomEntry(const float &_x, const float &_y, const float &_z, const int &_i)
  : x(_x), y(_y), z(_z), index(_i) {}
};

using namespace std;
extern "C" {
struct GridSearchPair {
  int ind1, ind2;
  GridSearchPair *next;
};

GridSearchPair *vmd_gridsearch3(const float *posA, int natomsA, const int *A,
                                const float *posB, int natomsB, const int *B,
                                float pairdist, int allow_double_counting, int maxpairs);

#define VMD_RAND_MAX 2147483647L
long vmd_random(void) {
#ifdef _MSC_VER
  return rand();
#else
  return random();
#endif
}

void vmd_srandom(unsigned int seed) {
#ifdef _MSC_VER
  srand(seed);
#else
  srandom(seed);
#endif
}

float distance2(const float *a, const float *b) {
  float delta = a[0] - b[0];
  float r2 = delta*delta;
  delta = a[1] - b[1];
  r2 += delta*delta;
  delta = a[2] - b[2];
  return r2 + delta*delta;
}

void vec_sub(float *a, const float *b, const float *c) {
  a[0]=b[0]-c[0];
  a[1]=b[1]-c[1];
  a[2]=b[2]-c[2];
}

void find_minmax(const float *pos, int n, const int *on,
                 float *min, float *max, int *oncount) {
  float x1, x2, y1, y2, z1, z2;
  int i, numon;

  // return immediately if there are no atoms, or no atoms are on.
  if (n < 1) return;

  // init on count, i = 0 because all atoms are on
  numon = 1;
  i = 0;
  // find first on atom EDIT: not needed
  // for (i=0; i<n; i++) {
  //   if (on[i]) {
  //     printf("on atom found\n");
  //     numon++;
  //     break;
  //   }
  // }

  if (i==n) {
    if (oncount != NULL)
      *oncount = numon;
    return;
  }

  // printf("x1 %f\n",pos[0]);
  // initialize min/max to first 'on' atom, and advance the counter.
  pos += 3*i;
  x1 = x2 = pos[0];
  y1 = y2 = pos[1];
  z1 = z2 = pos[2];
  pos += 3;
  i++;

  for (; i < n; i++) {
    // printf("x1: %f\n", x1);
    // if (on[i]) {
      if (pos[0] < x1) x1 = pos[0];
      else if (pos[0] > x2) x2 = pos[0];
      if (pos[1] < y1) y1 = pos[1];
      else if (pos[1] > y2) y2 = pos[1];
      if (pos[2] < z1) z1 = pos[2];
      else if (pos[2] > z2) z2 = pos[2];
      numon++;
    // }
    pos += 3;
  }
  min[0] = x1; min[1] = y1; min[2] = z1;
  max[0] = x2; max[1] = y2; max[2] = z2;

  if (oncount != NULL)
    *oncount = numon;
}


int find_minmax_selected(int n, const int *flgs, const float *pos,
                         float &_xmin, float &_ymin, float &_zmin,
                         float &_xmax, float &_ymax, float &_zmax) {
  int i;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  for (i=0; i<n; i++) if (flgs[i]) break;
  if (i==n) return FALSE;
  pos += 3*i;
  xmin=xmax=pos[0];
  ymin=ymax=pos[1];
  zmin=zmax=pos[2];
  pos += 3;
  i += 1;
  for (; i<n; i++, pos += 3) {
    if (!flgs[i]) continue;
    const float xi=pos[0];
    const float yi=pos[1];
    const float zi=pos[2];
    if (xmin>xi) xmin=xi;
    if (ymin>yi) ymin=yi;
    if (zmin>zi) zmin=zi;
    if (xmax<xi) xmax=xi;
    if (ymax<yi) ymax=yi;
    if (zmax<zi) zmax=zi;
  }
  _xmin=xmin;
  _ymin=ymin;
  _zmin=zmin;
  _xmax=xmax;
  _ymax=ymax;
  _zmax=zmax;
  return TRUE;
}



static void add_link(GridSearchPair *link, int i, int j) {
  link->next = (GridSearchPair *) malloc(sizeof(GridSearchPair));
  link->next->ind1 = i;
  link->next->ind2 = j;
  link->next->next = NULL;
}

int make_neighborlist(int **nbrlist, int xb, int yb, int zb) {
  int xi, yi, zi, aindex, xytotb;

  if (nbrlist == NULL)
    return -1;

  xytotb = xb * yb;
  aindex = 0;
  for (zi=0; zi<zb; zi++) {
    for (yi=0; yi<yb; yi++) {
      for (xi=0; xi<xb; xi++) {
        int nbrs[15]; // 14 neighbors, and a -1 to mark end of the list
        int n=0;
        nbrs[n++] = aindex;
        if (xi < xb-1) nbrs[n++] = aindex + 1;
        if (yi < yb-1) nbrs[n++] = aindex + xb;
        if (zi < zb-1) nbrs[n++] = aindex + xytotb;
        if (xi < (xb-1) && yi < (yb-1)) nbrs[n++] = aindex + xb + 1;
        if (xi < (xb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + 1;
        if (yi < (yb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + xb;
        if (xi < (xb-1) && yi > 0)      nbrs[n++] = aindex - xb + 1;
        if (xi > 0 && zi < (zb-1))     nbrs[n++] = aindex + xytotb - 1;
        if (yi > 0 && zi < (zb-1))     nbrs[n++] = aindex + xytotb - xb;
        if (xi < (xb-1) && yi < (yb-1) && zi < (zb-1))
                                       nbrs[n++] = aindex + xytotb + xb + 1;
        if (xi > 0 && yi < (yb-1) && zi < (zb-1))
                                       nbrs[n++] = aindex + xytotb + xb - 1;
        if (xi < (xb-1) && yi > 0 && zi < (zb-1))
                                       nbrs[n++] = aindex + xytotb - xb + 1;
        if (xi > 0 && yi > 0 && zi < (zb-1))
                                       nbrs[n++] = aindex + xytotb - xb - 1;
        nbrs[n++] = -1; // mark end of list

        int *lst = (int *) malloc(n*sizeof(int));
        if (lst == NULL)
          return -1; // return on failed allocations
        memcpy(lst, nbrs, n*sizeof(int));
        nbrlist[aindex] = lst;
        aindex++;
      }
    }
  }

  return 0;
}

GridSearchPair *vmd_gridsearch1(const float *pos,int natoms, const int *on,
                               float pairdist, int allow_double_counting, int maxpairs) {
  float min[3]={0,0,0}, max[3]={0,0,0};
  float sqdist;
  int i, j, xb, yb, zb, xytotb, totb, aindex;
  int **boxatom, *numinbox, *maxinbox, **nbrlist;
  int numon = 0;
  float sidelen[3], volume;
  int paircount = 0;
  int maxpairsreached = 0;
  sqdist = pairdist * pairdist;
  // printf("nc: %d \n",natoms);
  // find bounding box for selected atoms, and number of atoms in selection.
  find_minmax(pos, natoms, on, min, max, &numon);
  // printf("minmaxfound: %d,%f, min %f max %f\n",numon,sqdist,min[0],max[0]);
  // do sanity checks and complain if we've got bogus atom coordinates,
  // we shouldn't ever have density higher than 0.1 atom/A^3, but we'll
  // be generous and allow much higher densities.
  if (maxpairs != -1) {
    vec_sub(sidelen, max, min);
    // include estimate for atom radius (1 Angstrom) in volume determination
    volume = fabsf((sidelen[0] + 2.0f) * (sidelen[1] + 2.0f) * (sidelen[2] + 2.0f));
    if ((numon / volume) > 1.0) {
      // msgWarn << "vmd_gridsearch1: insane atom density" << sendmsg;
    }
  }

  // I don't want the grid to get too large, otherwise I could run out
  // of memory.  Octrees would be cool, but I'll just limit the grid size
  // and let the performance degrade a little for pathological systems.
  // Note that sqdist is what gets used for the actual distance checks;
  // from here on out pairdist is only used to set the grid size, so we
  // can set it to anything larger than the original pairdist.
  const int MAXBOXES = 4000000;
  totb = MAXBOXES + 1;

  float newpairdist = pairdist;
  float xrange = max[0]-min[0];
  float yrange = max[1]-min[1];
  float zrange = max[2]-min[2];
  do {
    // printf("pairdist: %f\n", pairdist);
    pairdist = newpairdist;
    const float invpairdist = 1.0f / pairdist;
    xb = ((int)(xrange*invpairdist))+1;
    yb = ((int)(yrange*invpairdist))+1;
    zb = ((int)(zrange*invpairdist))+1;
    xytotb = yb * xb;
    totb = xytotb * zb;
    newpairdist = pairdist * 1.26f; // cbrt(2) is about 1.26
  } while (totb > MAXBOXES || totb < 1); // check for integer wraparound too
  // printf("boxbuild\n");
  // printf("totb %d\n", totb);
  // 2. Sort each atom into appropriate bins
  boxatom = (int **) calloc(1, totb*sizeof(int *));
  numinbox = (int *) calloc(1, totb*sizeof(int));
  maxinbox = (int *) calloc(1, totb*sizeof(int));
  if (boxatom == NULL || numinbox == NULL || maxinbox == NULL) {
    if (boxatom != NULL)
      free(boxatom);
    if (numinbox != NULL)
      free(numinbox);
    if (maxinbox != NULL)
      free(maxinbox);
    // msgErr << "Gridsearch memory allocation failed, bailing out" << sendmsg;
    return NULL; // ran out of memory, bail out!
  }

  const float invpairdist = 1.0f / pairdist;
  for (i=0; i<natoms; i++) {
    // if (on[i]) {
      int axb, ayb, azb, aindex, num;

      // compute box index for new atom
      const float *loc = pos + 3*i;
      axb = (int)((loc[0] - min[0])*invpairdist);
      ayb = (int)((loc[1] - min[1])*invpairdist);
      azb = (int)((loc[2] - min[2])*invpairdist);

      // clamp box indices to valid range in case of FP error
      if (axb >= xb) axb = xb-1;
      if (ayb >= yb) ayb = yb-1;
      if (azb >= zb) azb = zb-1;

      aindex = azb * xytotb + ayb * xb + axb;

      // grow box if necessary
      if ((num = numinbox[aindex]) == maxinbox[aindex]) {
        boxatom[aindex] = (int *) realloc(boxatom[aindex], (num+4)*sizeof(int));
        maxinbox[aindex] += 4;
      }

      // store atom index in box
      boxatom[aindex][num] = i;
      numinbox[aindex]++;
    // }
  }
  free(maxinbox);

  nbrlist = (int **) calloc(1, totb*sizeof(int *));
  if (make_neighborlist(nbrlist, xb, yb, zb)) {
    if (boxatom != NULL) {
      for (i=0; i<totb; i++) {
        if (boxatom[i] != NULL) free(boxatom[i]);
      }
      free(boxatom);
    }
    if (nbrlist != NULL) {
      for (i=0; i<totb; i++) {
        if (nbrlist[i] != NULL) free(nbrlist[i]);
      }
      free(nbrlist);
    }
    free(numinbox);
    // msgErr << "Gridsearch memory allocation failed, bailing out" << sendmsg;
    return NULL; // ran out of memory, bail out!
  }

  // setup head of pairlist
  GridSearchPair *head, *cur;
  head = (GridSearchPair *) malloc(sizeof(GridSearchPair));
  head->next = NULL;
  cur = head;
  // printf("pairlist\n");
  // wkfmsgtimer *msgt = wkf_msg_timer_create(5);
  for (aindex = 0; (aindex < totb) && (!maxpairsreached); aindex++) {
    // printf("%d\n", aindex);
    int *tmpbox, *tmpnbr, *nbr;
    tmpbox = boxatom[aindex];
    tmpnbr = nbrlist[aindex];

    // if (wkf_msg_timer_timeout(msgt)) {
    //   char tmpbuf[128];
    //   sprintf(tmpbuf, "%6.2f", (100.0f * aindex) / (float) totb);
    //   msgInfo << "vmd_gridsearch1: " << tmpbuf << "% complete" << sendmsg;
    // }

    for (nbr = tmpnbr; (*nbr != -1) && (!maxpairsreached); nbr++) {
      int *nbrbox = boxatom[*nbr];
      for (i=0; (i<numinbox[aindex]) && (!maxpairsreached); i++) {
        int ind1 = tmpbox[i];
        // if (!on[ind1])
        //   continue;
        const float *p1 = pos + 3*ind1;
        int startj = 0;
        if (aindex == *nbr) startj = i+1;
        for (j=startj; (j<numinbox[*nbr]) && (!maxpairsreached); j++) {
          int ind2 = nbrbox[j];
          // if (on[ind2]) {
            const float *p2 = pos + 3*ind2;
            float ds2 = distance2(p1, p2);

            // ignore pairs between atoms with nearly identical coords
            if (ds2 < 0.001)
              continue;

            if (ds2 > sqdist)
              continue;

            if (maxpairs > 0) {
              if (paircount >= maxpairs) {
                maxpairsreached = 1;
                continue;
              }
            }

            add_link(cur, ind1, ind2);
            paircount++;

            // XXX double-counting still ignores atoms with same coords...
            if (allow_double_counting) {
              add_link(cur, ind2, ind1);
              paircount++;
            }
            cur = cur->next;
            cur->next = NULL;
          // }
        }
      }
    }
  }

  for (i=0; i<totb; i++) {
    free(boxatom[i]);
    free(nbrlist[i]);
  }
  free(boxatom);
  free(nbrlist);
  free(numinbox);

  cur = head->next;
  free(head);

  // if (maxpairsreached) {
  //   printf("maxpairs reached \n");
  // }

  return cur;
}

std::vector<std::vector<int>> find_contacts(const float *pos1, const float *pos2,
int nAtoms1, int nAtoms2, double cutoff) {
  std::vector<std::vector<int>> pairlist1(nAtoms1+nAtoms2);
  GridSearchPair *p, *tmp;
  int allowDouble = 1;
  int A[nAtoms1];
  fill_n(A, nAtoms1, 1);
  int B[nAtoms2];
  fill_n(B, nAtoms2, 1);
  GridSearchPair* pairs = vmd_gridsearch3(pos1, nAtoms1, A, pos2, nAtoms2, B, cutoff, allowDouble, -1);
  for (p = pairs; p != NULL; p = tmp) {
    int ind1=p->ind1;
    int ind2=p->ind2;
    pairlist1[ind1].push_back(ind2);
    // pairlist1[nAtoms1+ind2].push_back(ind1);
    tmp = p->next;
    free(p);
  }
  return pairlist1;
}


static double sasa_grid(const float *pos,int natoms, float pairdist, int allow_double_counting, int maxpairs, const float *radius,const int npts, double srad, int pointstyle, int restricted, const int* restrictedList) {
  int on[natoms];
  fill_n(on,natoms,1);
  // printf("natoms %d\n", natoms);
  // printf("pairdist %f\n", pairdist);
  // printf("maxpairs %d\n", maxpairs);
  GridSearchPair* pairs = vmd_gridsearch1(pos,natoms, on,  pairdist,  allow_double_counting,maxpairs);
  // vector<vector<int> > v(natoms,vector<int>());
  // printf("size: %d\n", v.size());
  ResizeArray<int> *pairlist = new ResizeArray<int>[natoms];
  GridSearchPair *p, *tmp;
  for (p = pairs; p != NULL; p = tmp) {
    int ind1=p->ind1;
    int ind2=p->ind2;
    // v[ind1].push_back(ind2);
    // v[ind2].push_back(ind1);
    pairlist[ind1].append(ind2);
    pairlist[ind2].append(ind1);
    tmp = p->next;
    free(p);
  }

  float *spherepts = new float[3*npts];

  if (pointstyle) {
    static const float RAND_MAX_INV = 1.0f/VMD_RAND_MAX;
    vmd_srandom(38572111);

    // All the spheres use the same random points.

    for (int i=0; i<npts; i++) {
      float u1 = (float) vmd_random();
      float u2 = (float) vmd_random();
      float z = 2.0f*u1*RAND_MAX_INV -1.0f;
      float phi = (float) (2.0f*M_PI*u2*RAND_MAX_INV);
      float R = sqrtf(1.0f-z*z);
      spherepts[3*i  ] = R*cosf(phi);
      spherepts[3*i+1] = R*sinf(phi);
      spherepts[3*i+2] = z;
      }
  }
  else {
    // All the spheres use the same spiral points.
    for (int k=0; k<npts; k++) {
      float h_k = 2.0 * (k - 1.0) / (npts - 1.0) - 1.0;
      float theta_k = acosf(h_k);
      float phi_k =0;
      if (! (k == 1 || k == npts)) {
          phi_k = fmod((phi_k + 3.6 / sqrt(npts * (1 - (h_k*h_k)))),(2 * M_PI));
      }
      spherepts[3*k  ] = cosf(phi_k) * sinf(theta_k);
      spherepts[3*k+1] = sinf(phi_k) * sinf(theta_k);
      spherepts[3*k+2] = cosf(theta_k);
    }
  }

  const float prefac = (float) (4 * M_PI / npts);
  float totarea = 0.0f;
  // compute area for each atom based on its pairlist
  for (int i = 0; i<natoms; i++) {
    // if (on[i]) {
      // only atoms in restrictsel contribute
      // printf("i: %d\n", restrictedList[i]);
      if (restricted && !restrictedList[i]) continue;
      const float *loc = pos+3*i;
      float rad = radius[i]+srad;
      float surfpos[3];
      int surfpts = npts;
      const ResizeArray<int> &nbrs = pairlist[i];
      // printf("neighbors: %d\n", nbrs.num());
      for (int j=0; j<npts; j++) {
        surfpos[0] = loc[0] + rad*spherepts[3*j  ];
        surfpos[1] = loc[1] + rad*spherepts[3*j+1];
        surfpos[2] = loc[2] + rad*spherepts[3*j+2];
        int ont = 1;
        for (int k=0; k<nbrs.num(); k++) {
          int ind = nbrs[k];
          const float *nbrloc = pos+3*ind;
          float radsq = radius[ind]+srad; radsq *= radsq;
          float dx = surfpos[0]-nbrloc[0];
          float dy = surfpos[1]-nbrloc[1];
          float dz = surfpos[2]-nbrloc[2];
          if (dx*dx + dy*dy + dz*dz <= radsq) {
            ont = 0;
            break;
          }
        }
        if (!ont) {
          surfpts--;
        }
      }
      float atomarea = prefac * rad * rad * surfpts;
      totarea += atomarea;
    // }
  }

  delete [] pairlist;
  delete [] spherepts;
  return totarea;
}



// PyObject* gridsearch(const float *pos,int natoms, float pairdist, int allow_double_counting, int maxpairs) {
//   int on[natoms];
//   fill_n(on,natoms,1);
//   printf("natoms %d\n", natoms);
//   printf("pairdist %f\n", pairdist);
//   printf("maxpairs %d\n", maxpairs);
//   GridSearchPair* pairs = vmd_gridsearch1(pos,natoms, on,  pairdist,  allow_double_counting,  maxpairs);
//   vector<vector<int> > v(natoms,vector<int>());
//   // printf("size: %d\n", v.size());
//   // ResizeArray<int> *pairlist = new ResizeArray<int>[natoms];
//   GridSearchPair *p, *tmp;
//   for (p = pairs; p != NULL; p = tmp) {
//     int ind1=p->ind1;
//     int ind2=p->ind2;
//     v[ind1].push_back(ind2);
//     v[ind2].push_back(ind1);
//     // printf("%f\n", v[ind2]);
//     // pairlist[ind1].append(ind2);
//     // pairlist[ind2].append(ind1);
//     // printf("%d\n", ind1);
//     tmp = p->next;
//     free(p);
//   }
//   PyGILState_STATE gstate = PyGILState_Ensure();
//   PyObject* result = PyList_New(0);
//   vector< vector<int> >::iterator row;
//   vector<int>::iterator col;
//   for (row = v.begin(); row != v.end(); row++) {
//     // printf("row %d\n", *row);
//     PyObject* tempo = PyList_New(0);
//     for (col = row->begin(); col != row->end(); col++) {
//         // printf("%d \n", *col);
//         PyList_Append(tempo, PyInt_FromLong(*col));
//     }
//     PyList_Append(result,tempo);
//     tempo = NULL;
//   }
//   PyGILState_Release(gstate);
//   return result;
// }

}

typedef ResizeArray<AtomEntry> atomlist;
struct FindWithinData {
  int nthreads;
  int tid;
  int totb;
  int xytotb;
  int xb;
  int yb;
  int zb;
  float r2;
  const float * xyz;
  const atomlist * flgatoms;
  const atomlist * otheratoms;
  int * flgs;
  FindWithinData() : flgatoms(0), otheratoms(0), flgs(0) {}
  ~FindWithinData() { if (flgs) free(flgs); }
};

#define MAXGRIDDIM 31
static int* find_within(const float *xyz, int *flgs, int *others, int num, float r) {
  int i;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  float oxmin, oymin, ozmin, oxmax, oymax, ozmax;
  float xwidth, ywidth, zwidth;
  const float *pos;

  int *result = new int[num];
  fill_n(result,num,0);

  // for (size_t i = 0; i < num; i++) {
  //   pos=xyz+3*i;
  //   printf("%f %f %f\n", pos[0],pos[1],pos[2]);
  // }

  // find min/max bounds of atom coordinates in flgs
  if (!find_minmax_selected(num, flgs, xyz, xmin, ymin, zmin, xmax, ymax, zmax) ||
    !find_minmax_selected(num, others, xyz, oxmin, oymin, ozmin, oxmax, oymax, ozmax)) {
    memset(flgs, 0, num*sizeof(int));
    return result;
  }

  // Find the set of atoms with the smallest extent; here we use the sum
  // of the box dimensions though other choices might be better.
  float size = xmax+ymax+zmax - (xmin+ymin+zmin);
  float osize = oxmax+oymax+ozmax - (oxmin+oymin+ozmin);
  if (osize < size) {
    xmin=oxmin;
    ymin=oymin;
    zmin=ozmin;
    xmax=oxmax;
    ymax=oymax;
    zmax=ozmax;
  }

  // Generate a grid of mesh size r based on the computed size of the molecule.
  // We limit the size of the grid cell dimensions so that we don't get too
  // many grid cells.
  xwidth = (xmax-xmin)/(MAXGRIDDIM-1);
  if (xwidth < r) xwidth = r;
  ywidth = (ymax-ymin)/(MAXGRIDDIM-1);
  if (ywidth < r) ywidth = r;
  zwidth = (zmax-zmin)/(MAXGRIDDIM-1);
  if (zwidth < r) zwidth = r;

  // Adjust the bounds so that we include atoms that are in the outermost
  // grid cells.
  xmin -= xwidth;
  xmax += xwidth;
  ymin -= ywidth;
  ymax += ywidth;
  zmin -= zwidth;
  zmax += zwidth;

  // The number of grid cells needed in each dimension is
  // (int)((xmax-xmin)/xwidth) + 1
  const int xb = (int)((xmax-xmin)/xwidth) + 1;
  const int yb = (int)((ymax-ymin)/ywidth) + 1;
  const int zb = (int)((zmax-zmin)/zwidth) + 1;

  int xytotb = yb * xb;
  int totb = xytotb * zb;

  atomlist* flgatoms   = new atomlist[totb];
  atomlist* otheratoms = new atomlist[totb];

  float ixwidth = 1.0f/xwidth;
  float iywidth = 1.0f/ywidth;
  float izwidth = 1.0f/zwidth;
  for (i=0; i<num; i++) {
    if (!flgs[i] && !others[i]) continue;
    pos=xyz+3*i;
    float xi = pos[0];
    float yi = pos[1];
    float zi = pos[2];
    if (xi<xmin || xi>xmax || yi<ymin || yi>ymax || zi<zmin || zi>zmax) {
      continue;
    }
    AtomEntry entry(xi,yi,zi,i);
    int axb = (int)((xi - xmin)*ixwidth);
    int ayb = (int)((yi - ymin)*iywidth);
    int azb = (int)((zi - zmin)*izwidth);

    // Due to floating point error in the calcuation of bin widths, we
    // have to range clamp the computed box indices.
    if (axb==xb) axb=xb-1;
    if (ayb==yb) ayb=yb-1;
    if (azb==zb) azb=zb-1;

    int aindex = azb*xytotb + ayb*xb + axb;
    // TODO maybe change so that others does not contain the atoms of flags
    if (others[i]) otheratoms[aindex].append(entry);
    if (  flgs[i])   flgatoms[aindex].append(entry);
  }

  memset(flgs, 0, num*sizeof(int));
  const float r2 = (float) (r*r);

  // set up workspace for multithreaded calculation
  int nthreads;
#ifdef VMDTHREADS
  nthreads = wkf_thread_numprocessors();
  wkf_thread_t * threads = (wkf_thread_t *)calloc(nthreads, sizeof(wkf_thread_t));
#else
  nthreads = 1;
#endif
  FindWithinData *data = new FindWithinData[nthreads];
  for (i=0; i<nthreads; i++) {
    data[i].nthreads = nthreads;
    data[i].tid = i;
    data[i].totb = totb;
    data[i].xytotb = xytotb;
    data[i].xb = xb;
    data[i].yb = yb;
    data[i].zb = zb;
    data[i].r2 = r2;
    data[i].xyz = xyz;
    data[i].flgatoms = flgatoms;
    data[i].otheratoms = otheratoms;
    data[i].flgs = (int *)calloc(num, sizeof(int));
  }
#ifdef VMDTHREADS
  for (i=0; i<nthreads; i++) {
    wkf_thread_create(threads+i, find_within_routine, data+i);
  }
  for (i=0; i<nthreads; i++) {
    wkf_thread_join(threads[i], NULL);
  }
  free(threads);
#else
  find_within_routine(data);
#endif

  // combine results
  for (i=0; i<nthreads; i++) {
    const int *tflg = data[i].flgs;
    for (int j=0; j<num; j++) {
      flgs[j] |= tflg[j];
      result[j] = tflg[j];
      // printf("j: %d , %d\n", j, flgs[j]);
    }
  }
  delete [] data;
  delete [] flgatoms;
  delete [] otheratoms;
  return result;
}

void * find_within_routine( void *v ) {
  FindWithinData *data = (FindWithinData *)v;
  const int nthreads = data->nthreads;
  const int tid = data->tid;
  const int totb = data->totb;
  const int xytotb = data->xytotb;
  const int xb = data->xb;
  const int yb = data->yb;
  const int zb = data->zb;
  const float r2 = data->r2;
  const atomlist * flgatoms = data->flgatoms;
  const atomlist * otheratoms = data->otheratoms;
  int * flgs = data->flgs;

  // Loop over boxes, checking for flg atoms and other atoms within one
  // box of each other.  When one is found, mark the flag.
  for (int aindex = tid; aindex<totb; aindex += nthreads) {
    // Figure out the neighbors for this box
    int zi = aindex/xytotb;
    int ytmp = aindex - zi*xytotb;
    int yi = ytmp/xb;
    int xi = ytmp - yi*xb;
    int nbrs[14];
    int n=0;
    nbrs[n++] = aindex;     // Always include self
    if (xi < xb-1) nbrs[n++] = aindex + 1;
    if (yi < yb-1) nbrs[n++] = aindex + xb;
    if (zi < zb-1) nbrs[n++] = aindex + xytotb;
    if (xi < (xb-1) && yi < (yb-1)) nbrs[n++] = aindex + xb + 1;
    if (xi < (xb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + 1;
    if (yi < (yb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + xb;
    if (xi < (xb-1) && yi > 0)      nbrs[n++] = aindex - xb + 1;
    if (xi > 0 && zi < (zb-1))     nbrs[n++] = aindex + xytotb - 1;
    if (yi > 0 && zi < (zb-1))     nbrs[n++] = aindex + xytotb - xb;
    if (xi < (xb-1) && yi < (yb-1) && zi < (zb-1))
                                   nbrs[n++] = aindex + xytotb + xb + 1;
    if (xi > 0 && yi < (yb-1) && zi < (zb-1))
                                   nbrs[n++] = aindex + xytotb + xb - 1;
    if (xi < (xb-1) && yi > 0 && zi < (zb-1))
                                   nbrs[n++] = aindex + xytotb - xb + 1;
    if (xi > 0 && yi > 0 && zi < (zb-1))
                                   nbrs[n++] = aindex + xytotb - xb - 1;

    const atomlist& boxflg = flgatoms[aindex];
    // Compare the atoms in boxflg to those in nbrother
    int i;
    for (i=0; i<boxflg.num(); i++) {
      const AtomEntry &flgentry = boxflg[i];
      int flgind = flgentry.index;

      // Compare flag atoms in this box to other atoms in neighbor boxes,
      for (int inbr=0; inbr<n; inbr++) {  // Loop over neighbors
        if (flgs[flgind]) break;

        // Fetch a neighboring otheratoms to compare to boxflg
        int nbr = nbrs[inbr];
        const atomlist& nbrother = otheratoms[nbr];
        for (int j=0; j<nbrother.num(); j++) {
          const AtomEntry &otherentry = nbrother[j];
          float dx2 = flgentry.x - otherentry.x; dx2 *= dx2;
          float dy2 = flgentry.y - otherentry.y; dy2 *= dy2;
          float dz2 = flgentry.z - otherentry.z; dz2 *= dz2;
          if (dx2 + dy2 + dz2 < r2) {
            flgs[flgind] = 1;
            break;
          }
        }
      }
    }

    // compare other atoms in this box to flag atoms in the neighbors.
    const atomlist& boxother = otheratoms[aindex];

    for (int inbr=0; inbr<n; inbr++) {  // Loop over neighbors
      int nbr = nbrs[inbr];

      // Fetch a neighboring flgatoms to compare to boxother
      const atomlist& nbrflg = flgatoms[nbr];

      // Compare the atoms in boxother to those in nbrflg
      for (i=0; i<nbrflg.num(); i++) {
        const AtomEntry &flgentry = nbrflg[i];
        int flgind = flgentry.index;
        // The next test helps a lot when the boxes are large, but hurts
        // a little when the boxes are small.
        if (flgs[flgind]) continue;
        for (int j=0; j<boxother.num(); j++) {
          const AtomEntry &otherentry = boxother[j];
          float dx2 = flgentry.x - otherentry.x; dx2 *= dx2;
          float dy2 = flgentry.y - otherentry.y; dy2 *= dy2;
          float dz2 = flgentry.z - otherentry.z; dz2 *= dz2;
          if (dx2 + dy2 + dz2 < r2) {
            flgs[flgind] = 1;
            break;
          }
        }
      }
    }
  }
  return NULL;
}

static int make_neighborlist_sym(int **nbrlist, int xb, int yb, int zb) {
  int xi, yi, zi, aindex, xytotb;

  if (nbrlist == NULL)
    return -1;

  xytotb = xb * yb;
  aindex = 0;
  for (zi=0; zi<zb; zi++) {
    for (yi=0; yi<yb; yi++) {
      for (xi=0; xi<xb; xi++) {
        int nbrs[28]; // 27 neighbors, and a -1 to mark end of the list
        int n=0;
        nbrs[n++] = aindex;
        if (xi < xb-1) nbrs[n++] = aindex + 1;
        if (yi < yb-1) nbrs[n++] = aindex + xb;
        if (zi < zb-1) nbrs[n++] = aindex + xytotb;
        if (xi < (xb-1) && yi < (yb-1)) nbrs[n++] = aindex + xb + 1;
        if (xi < (xb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + 1;
        if (yi < (yb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + xb;
        if (xi < (xb-1) && yi) nbrs[n++] = aindex - xb + 1;
        if (xi && zi < (zb-1)) nbrs[n++] = aindex + xytotb - 1;
        if (yi && zi < (zb-1)) nbrs[n++] = aindex + xytotb - xb;
        if (xi < (xb-1) && yi < (yb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + xb + 1;
        if (xi && yi < (yb-1) && zi < (zb-1)) nbrs[n++] = aindex + xytotb + xb - 1;
        if (xi < (xb-1) && yi && zi < (zb-1)) nbrs[n++] = aindex + xytotb - xb + 1;
        if (xi && yi && zi < (zb-1)) nbrs[n++] = aindex + xytotb - xb - 1;

        if (xi) nbrs[n++] = aindex - 1;
        if (yi) nbrs[n++] = aindex - xb;
        if (zi) nbrs[n++] = aindex - xytotb;
        if (xi && yi) nbrs[n++] = aindex - xb - 1;
        if (xi && zi) nbrs[n++] = aindex - xytotb - 1;
        if (yi && zi) nbrs[n++] = aindex - xytotb - xb;
        if (xi && yi < (yb-1)) nbrs[n++] = aindex + xb - 1;
        if (xi < (xb-1) && zi) nbrs[n++] = aindex - xytotb + 1;
        if (yi < (yb-1) && zi) nbrs[n++] = aindex - xytotb + xb;
        if (xi && yi && zi) nbrs[n++] = aindex - xytotb - xb - 1;
        if (xi < (xb-1) && yi && zi) nbrs[n++] = aindex - xytotb - xb + 1;
        if (xi && yi < (yb-1) && zi) nbrs[n++] = aindex - xytotb + xb - 1;
        if (xi < (xb-1) && yi < (yb-1) && zi) nbrs[n++] = aindex - xytotb + xb + 1;
        nbrs[n++] = -1; // mark end of list

        int *lst = (int *) malloc(n*sizeof(int));
        if (lst == NULL)
          return -1; // return on failed allocations
        memcpy(lst, nbrs, n*sizeof(int));
        nbrlist[aindex] = lst;
        aindex++;
      }
    }
  }

  return 0;
}


GridSearchPair *vmd_gridsearch3(const float *posA, int natomsA, const int *A,
                                const float *posB, int natomsB, const int *B,
                                float pairdist, int allow_double_counting, int maxpairs) {

  if (!natomsA || !natomsB) return NULL;

  if (allow_double_counting == -1) { //default
    if (posA == posB && natomsA == natomsB)
      allow_double_counting = FALSE;
    else
      allow_double_counting = TRUE;
  }

  // if same molecule and *A[] == *B[], it is a lot faster to use gridsearch1
  if (posA == posB && natomsA == natomsB) {
    bool is_equal = TRUE;
    for (int i=0; i<natomsA && is_equal; i++) {
      if (A[i] != B[i])
        is_equal = FALSE;
    }
    if (is_equal)
      return vmd_gridsearch1(posA, natomsA, A, pairdist, allow_double_counting, maxpairs);
  }

  float min[3], max[3], sqdist;
  float minB[3], maxB[3]; //tmp storage
  int i, j, xb, yb, zb, xytotb, totb, aindex;
  int **boxatomA, *numinboxA, *maxinboxA;
  int **boxatomB, *numinboxB, *maxinboxB;
  int **nbrlist;
  float sidelen[3], volume;
  int numonA = 0;   int numonB = 0;
  int paircount = 0;
  int maxpairsreached = 0;
  sqdist = pairdist * pairdist;

  // 1. Find grid size for binning
  // find bounding box for selected atoms
  find_minmax(posA, natomsA, A, min, max, &numonA);
  find_minmax(posB, natomsB, B, minB, maxB, &numonB);

  // If no atoms were selected we don't have to go on
  if (!numonA || !numonB) return NULL;

  for (i=0; i<3; i++) {
    if (minB[i] < min[i]) min[i] = minB[i];
    if (maxB[i] > max[i]) max[i] = maxB[i];
  }

  // do sanity checks and complain if we've got bogus atom coordinates,
  // we shouldn't ever have density higher than 0.1 atom/A^3, but we'll
  // be generous and allow much higher densities.
  if (maxpairs != -1) {
    vec_sub(sidelen, max, min);
    // include estimate for atom radius (1 Angstrom) in volume determination
    volume = fabsf((sidelen[0] + 2.0f) * (sidelen[1] + 2.0f) * (sidelen[2] + 2.0f));
    if (((numonA + numonB) / volume) > 1.0) {
      // msgWarn << "vmd_gridsearch3: insane atom density" << sendmsg;
    }
  }

  // I don't want the grid to get too large, otherwise I could run out
  // of memory.  Octrees would be cool, but I'll just limit the grid size
  // and let the performance degrade a little for pathological systems.
  // Note that sqdist is what gets used for the actual distance checks;
  // from here on out pairdist is only used to set the grid size, so we
  // can set it to anything larger than the original pairdist.
  const int MAXBOXES = 4000000;
  totb = MAXBOXES + 1;

  float newpairdist = pairdist;
  float xrange = max[0]-min[0];
  float yrange = max[1]-min[1];
  float zrange = max[2]-min[2];
  do {
    pairdist = newpairdist;
    const float invpairdist = 1.0f / pairdist;
    xb = ((int)(xrange*invpairdist))+1;
    yb = ((int)(yrange*invpairdist))+1;
    zb = ((int)(zrange*invpairdist))+1;
    xytotb = yb * xb;
    totb = xytotb * zb;
    newpairdist = pairdist * 1.26f; // cbrt(2) is about 1.26
  } while (totb > MAXBOXES || totb < 1); // check for integer wraparound too

  // 2. Sort each atom into appropriate bins
  boxatomA = (int **) calloc(1, totb*sizeof(int *));
  numinboxA = (int *) calloc(1, totb*sizeof(int));
  maxinboxA = (int *) calloc(1, totb*sizeof(int));
  if (boxatomA == NULL || numinboxA == NULL || maxinboxA == NULL) {
    if (boxatomA != NULL)
      free(boxatomA);
    if (numinboxA != NULL)
      free(numinboxA);
    if (maxinboxA != NULL)
      free(maxinboxA);
    // msgErr << "Gridsearch memory allocation failed, bailing out" << sendmsg;
    return NULL; // ran out of memory, bail out!
  }

  const float invpairdist = 1.0f / pairdist;
  for (i=0; i<natomsA; i++) {
    if (A[i]) {
      int axb, ayb, azb, aindex, num;

      // compute box index for new atom
      const float *loc = posA + 3L*i;
      axb = (int)((loc[0] - min[0])*invpairdist);
      ayb = (int)((loc[1] - min[1])*invpairdist);
      azb = (int)((loc[2] - min[2])*invpairdist);

      // clamp box indices to valid range in case of FP error
      if (axb >= xb) axb = xb-1;
      if (ayb >= yb) ayb = yb-1;
      if (azb >= zb) azb = zb-1;

      aindex = azb * xytotb + ayb * xb + axb;

      // grow box if necessary
      if ((num = numinboxA[aindex]) == maxinboxA[aindex]) {
        boxatomA[aindex] = (int *) realloc(boxatomA[aindex], (num+4)*sizeof(int));
        maxinboxA[aindex] += 4;
      }

      // store atom index in box
      boxatomA[aindex][num] = i;
      numinboxA[aindex]++;
    }
  }
  free(maxinboxA);

  boxatomB = (int **) calloc(1, totb*sizeof(int *));
  numinboxB = (int *) calloc(1, totb*sizeof(int));
  maxinboxB = (int *) calloc(1, totb*sizeof(int));
  if (boxatomB == NULL || numinboxB == NULL || maxinboxB == NULL) {
    if (boxatomB != NULL)
      free(boxatomB);
    if (numinboxB != NULL)
      free(numinboxB);
    if (maxinboxB != NULL)
      free(maxinboxB);
    // msgErr << "Gridsearch memory allocation failed, bailing out" << sendmsg;
    return NULL; // ran out of memory, bail out!
  }

  for (i=0; i<natomsB; i++) {
    if (B[i]) {
      int axb, ayb, azb, aindex, num;

      // compute box index for new atom
      const float *loc = posB + 3L*i;
      axb = (int)((loc[0] - min[0])*invpairdist);
      ayb = (int)((loc[1] - min[1])*invpairdist);
      azb = (int)((loc[2] - min[2])*invpairdist);

      // clamp box indices to valid range in case of FP error
      if (axb >= xb) axb = xb-1;
      if (ayb >= yb) ayb = yb-1;
      if (azb >= zb) azb = zb-1;

      aindex = azb * xytotb + ayb * xb + axb;

      // grow box if necessary
      if ((num = numinboxB[aindex]) == maxinboxB[aindex]) {
        boxatomB[aindex] = (int *) realloc(boxatomB[aindex], (num+4)*sizeof(int));
        maxinboxB[aindex] += 4;
      }

      // store atom index in box
      boxatomB[aindex][num] = i;
      numinboxB[aindex]++;
    }
  }
  free(maxinboxB);


  // 3. Build pairlists of atoms less than sqrtdist apart
  nbrlist = (int **) calloc(1, totb*sizeof(int *));
  if (make_neighborlist_sym(nbrlist, xb, yb, zb)) {
    if (boxatomA != NULL) {
      for (i=0; i<totb; i++) {
        if (boxatomA[i] != NULL) free(boxatomA[i]);
      }
      free(boxatomA);
    }
    if (boxatomB != NULL) {
      for (i=0; i<totb; i++) {
        if (boxatomB[i] != NULL) free(boxatomB[i]);
      }
      free(boxatomB);
    }
    if (nbrlist != NULL) {
      for (i=0; i<totb; i++) {
        if (nbrlist[i] != NULL) free(nbrlist[i]);
      }
      free(nbrlist);
    }
    free(numinboxA);
    free(numinboxB);
    // msgErr << "Gridsearch memory allocation failed, bailing out" << sendmsg;
    return NULL; // ran out of memory, bail out!
  }

  // setup head of pairlist
  GridSearchPair *head, *cur;
  head = (GridSearchPair *) malloc(sizeof(GridSearchPair));
  head->next = NULL;
  cur = head;

  for (aindex = 0; aindex < totb; aindex++) {
    int *tmpbox, *tmpnbr, *nbr;
    tmpbox = boxatomA[aindex];
    tmpnbr = nbrlist[aindex];
    for (nbr = tmpnbr; (*nbr != -1) && (!maxpairsreached); nbr++) {
      int *nbrboxB = boxatomB[*nbr];
      for (i=0; (i<numinboxA[aindex]) && (!maxpairsreached); i++) {
        const float *p1;
        int ind1 = tmpbox[i];
        p1 = posA + 3L*ind1;
        for (j=0; (j<numinboxB[*nbr]) && (!maxpairsreached); j++) {
          const float *p2;
          int ind2 = nbrboxB[j];
          p2 = posB + 3L*ind2;
          if (!allow_double_counting && B[ind1] && A[ind2] && ind2<=ind1) continue; //don't double-count bonds XXX
          if (distance2(p1,p2) > sqdist) continue;

          if (maxpairs > 0) {
            if (paircount >= maxpairs) {
              maxpairsreached = 1;
              continue;
            }
          }

          add_link(cur, ind1, ind2);
          paircount++;
          cur = cur->next;
          cur->next = NULL;
        }
      }
    }
  }

  for (i=0; i<totb; i++) {
    free(boxatomA[i]);
    free(boxatomB[i]);
    free(nbrlist[i]);
  }
  free(boxatomA);
  free(boxatomB);
  free(nbrlist);
  free(numinboxA);
  free(numinboxB);

  cur = head->next;
  free(head);

  // if (maxpairsreached)
    // msgErr << "gridsearch3: exceeded pairlist sanity check, aborted" << sendmsg;

  return cur;
}
