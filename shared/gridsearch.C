#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ResizeArray.h"
#include <vector>
// #include "utilities.h"

using namespace std;
extern "C" {
struct GridSearchPair {
  int ind1, ind2;
  GridSearchPair *next;
};

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

  printf("x1 %f\n",pos[0]);
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
  printf("nc: %d \n",natoms);
  // find bounding box for selected atoms, and number of atoms in selection.
  find_minmax(pos, natoms, on, min, max, &numon);
  printf("minmaxfound: %d,%f, min %f max %f\n",numon,sqdist,min[0],max[0]);
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
  printf("boxbuild\n");
  printf("totb %d\n", totb);
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
  printf("pairlist\n");
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


double sasa_grid(const float *pos,int natoms, float pairdist, int allow_double_counting, int maxpairs, const float *radius,const int npts, double srad, int pointstyle, int restricted, const int* restrictedList) {
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
