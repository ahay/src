/* Applying tapering absorbing boundary condition for 2d wavefield (stored as oned array) */
/*
  Copyright (C) 2014 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "absorb.h"

static int nx, nz, nx2, nz2, nbt, nbb, nbl, nbr;
static float ct, cb, cl, cr;
static float *wt, *wb, *wl, *wr;

void abc_init(int n1,  int n2    /*model size*/,
	      int n12, int n22   /*padded model size*/,
	      int nb1, int nb2   /*top, bottom*/,
	      int nb3, int nb4   /*left, right*/,
	      float c1, float c2 /*top, bottom*/,
	      float c3, float c4 /*left, right*/)
/*< initialization >*/
{
    int c;
    nz = n1;
    nx = n2;
    nz2= n12;
    nx2= n22;
    nbt = nb1;
    nbb = nb2;
    nbl = nb3;
    nbr = nb4;
    ct = c1;
    cb = c2;
    cl = c3;
    cr = c4;
    if(nbt) wt =  sf_floatalloc(nbt);
    if(nbb) wb =  sf_floatalloc(nbb);
    if(nbl) wl =  sf_floatalloc(nbl);
    if(nbr) wr =  sf_floatalloc(nbr);
    c=0;
    abc_cal(c,nbt,ct,wt);
    abc_cal(c,nbb,cb,wb);
    abc_cal(c,nbl,cl,wl);
    abc_cal(c,nbr,cr,wr);
}
   

void abc_close(void)
/*< free memory allocation>*/
{
    if(nbt) free(wt);
    if(nbb) free(wb);
    if(nbl) free(wl);
    if(nbr) free(wr);
}

void abc_apply(float *a /*2-D matrix*/) 
/*< boundary decay>*/
{
    int i;
    int iz, ix;

    /* top */
#ifdef _OPENMP
#pragma omp parallel default(shared) private(iz,ix,i)
{
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbt; iz++) {  
        for (ix=nbl; ix < nx-nbr; ix++) {
	  i = nz2*ix + iz;
	  a[i] *= wt[iz];
        }
    }
    /* bottom */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbb; iz++) {  
        for (ix=nbl; ix < nx-nbr; ix++) {
	  i = nz2*ix + nz-1-iz;
	  a[i] *= wb[iz];
        }
    }
    /* left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=nbt; iz < nz-nbb; iz++) {  
        for (ix=0; ix < nbl; ix++) {
	  i = nz2*ix + iz;
	  a[i] *= wl[ix];
        }
    }
    /* right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=nbt; iz < nz-nbb; iz++) {  
        for (ix=0; ix < nbr; ix++) {
	  i = nz2*(nx-1-ix) + iz;
          a[i] *= wr[ix];
        }
    }
    /* top left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbt; iz++) {  
        for (ix=0; ix < nbl; ix++) {
	  i = nz2*ix + iz;
	  a[i] *= iz>ix? wl[ix]:wt[iz];
        }
    }
    /* top right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbt; iz++) {  
        for (ix=0; ix < nbr; ix++) {
	  i = nz2*(nx-1-ix) + iz;
	  a[i] *= iz>ix? wr[ix]:wt[iz];
        }
    }
    /* bottom left */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbb; iz++) {  
        for (ix=0; ix < nbl; ix++) {
	  i = nz2*ix + nz-1-iz;
          a[i] *= iz>ix? wl[ix]:wb[iz];
        }
    }
    /* bottom right */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nbb; iz++) {  
        for (ix=0; ix < nbr; ix++) {
	  i = nz2*(nx-1-ix) + nz-1-iz;
          a[i] *= iz>ix? wr[ix]:wb[iz];
        }
    }

#ifdef _OPENMP
}
#endif
}

void abc_cal(int abc /* decaying type*/,
             int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< find absorbing coefficients >*/
{
    int ib;
    /*const float pi=SF_PI;*/
    if(!nb) return;
    switch(abc) {
    default:
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ib)
#endif
        for(ib=0; ib<nb; ib++){
	    w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
	}
    }
}
