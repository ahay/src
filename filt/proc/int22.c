#include <math.h>
#include "int22.h"

#include "bool.h"
#include "adjnull.h"
#include "seperr.h"


#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif


static int nd, nf, m[2];
static int *nxy[2];
static bool *mask;
static float **w1, **w2;

void  int22_init (float** coord, float o[], float d[], int n[], 
		  interpolator interp, int nf_in, int nd_in)
{
  int id, i; 
  float x[2], rx;

  nf = nf_in;
  nd = nd_in;
  m[0] = n[0];
  m[1] = n[1];

  nxy[0] = (int*) sepalloc(nd,sizeof(int));
  nxy[1] = (int*) sepalloc(nd,sizeof(int));
  mask = (bool*) sepalloc(nd,sizeof(bool));
  w1 = (float**) sepalloc2(nf,nd,sizeof(float));
  w2 = (float**) sepalloc2(nf,nd,sizeof(float));

  for (id = 0; id < nd; id++) {
    for (i =0; i < 2; i++) {
      rx = (coord[i][id] - o[i])/d[i];
      nxy[i][id] = (int) floor(rx + 1. - 0.5*nf);
      x[i] = rx - floor(rx);
    }
    
    if (nxy[0][id] > - nf && nxy[0][id] < n[0] &&
	nxy[1][id] > - nf && nxy[1][id] < n[1]) {
      mask[id] = FALSE; 
      interp (x[0], nf, w1[id]);
      interp (x[1], nf, w2[id]);
    } else {
      mask[id] = TRUE;
      for (i =0; i < nf; i++) {
	w1[id][i] = 0.; 
	w2[id][i] = 0.;
      } 
    }
  }
}

void  int22_lop (bool adj, bool add, int nm, int nd, float* x, float* ord)
{ 
  int id, i0, j0, i, j, im;
  float w;

  adjnull (adj,add,nm,nd,x,ord);

  for (id=0; id < nd; id++) {
    if (mask[id] == TRUE) continue;
    i0 = nxy[0][id]; 
    j0 = nxy[1][id]; 
    for (j = MAX(0,-j0); j < MIN(nf,m[1]-j0); j++) {
      w = w2[id][j];
      for (i = MAX(0,-i0); i < MIN(nf,m[0]-i0); i++) { 
	im = (i+i0) + (j+j0)*m[0];
	if( adj) { 
	  x[im] += ord[id] * w * w1[id][i];
	} else {
	  ord[id] += x[im] * w * w1[id][i];
	}
      }
    }
  }
}

void int22_close (void)
{
  free (nxy[0]);
  free (nxy[1]);
  free (mask);
  free2 (nd,(void**)w1);
  free2 (nd,(void**)w2);
}
