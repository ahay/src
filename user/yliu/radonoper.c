/* Radon operator. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include <math.h>
#include <rsf.h>
#include "radonoper.h"
#include "radon.h"     

static int nt;                  /* number of time samples */
static int nt2;                 /* extended number of time samples */
static int nx;                  /* number of offsets */ 
static int np;                  /* number of slopes */
static int nw;                  /* number of frequencies */
static float dw;                /* frequency increment */

static kiss_fftr_cfg forw, invs;
static sf_complex *dd;           /* data (CMP gather)    */
static sf_complex *mm;           /* model (Radon gather) */
static sf_complex **cm, **cd;    /* model and data storage */
static float *tt;                /* trace */
static float *xx;                /* offset header */

void radonoper_init (int nt_in, float dt_in, float t0_in,
		     int nx_in, float dx_in, float ox_in, float x0_in,
		     int np_in, float dp_in, float p0_in,
		     bool par_in)
/*< initialize >*/
{
    int ix;
    nx = nx_in;
    nt = nt_in;
    np = np_in;
    nt2 = 2*kiss_fft_next_fast_size(nt_in);
    nw  = nt2/2+1;
    dw  = 2.0*SF_PI/(nt2*dt_in);

    dd = sf_complexalloc(nx_in);
    mm = sf_complexalloc(np_in);  
    tt = sf_floatalloc (nt2);
    cm = sf_complexalloc2 (nw,np_in);
    cd = sf_complexalloc2 (nw,nx_in);   

    forw = kiss_fftr_alloc(nt2,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt2,1,NULL,NULL);
    dd = sf_complexalloc(nx_in);
    mm = sf_complexalloc(np_in);   
    xx = sf_floatalloc(nx_in);

    for (ix=0; ix < nx_in; ix++) {
	xx[ix] = ox_in + ix*dx_in;
    }
    for (ix=0; ix < nx_in; ix++) {
	if (par_in) {
	    xx[ix] *= xx[ix]/(x0_in*x0_in); 
	} else if (1. != x0_in) {
	    xx[ix] /= x0_in;
	}
    }
 
    radon_init (nx_in, np_in, dp_in, p0_in); /* initiliaze radon operator */
    
}

void radonoper_lop (bool adj, bool add, int nxx, int nyy, float *x, float *y)
/*< linear operator >*/
{
    int ix, it, ip, iw;
    float w;

    sf_adjnull (adj,add,nxx,nyy,x,y);

    if (adj) {
	for (ix=0; ix < nx; ix++) { /* loop over offsets */
	    for (it=0; it < nt; it++) {
		tt[it] = y[ix*nt+it];
	    }
	    for (it=nt; it < nt2; it++) {
		tt[it] = 0.;
	    }
	    /* FFT to frequency */
	    kiss_fftr(forw,tt, (kiss_fft_cpx *) cd[ix]);
	}
    } else { /* modeling */
	for (ip=0; ip < np; ip++) { /* loop over slopes */
	    for (it=0; it < nt; it++) {
		tt[it] = x[ip*nt+it];
	    }
	    for (it=nt; it < nt2; it++) {
		tt[it] = 0.;
	    }

	    /* FFT to frequency */
	    kiss_fftr(forw,tt, (kiss_fft_cpx *) cm[ip]);
	}
    }
	
    for (iw=0; iw < nw; iw++) { /* loop over frequencies */
	w = iw*dw;
	
	if (adj) {
	    for (ix=0; ix < nx; ix++) { /* loop over offsets */
#ifdef SF_HAS_COMPLEX_H
		dd[ix] = cd[ix][iw]/nt2; /* transpose */
#else
		dd[ix] = sf_crmul(cd[ix][iw],1.0/nt2); 
#endif
	    } 
	} else {
	    for (ip=0; ip < np; ip++) { /* loop over slopes */
#ifdef SF_HAS_COMPLEX_H
		mm[ip] = cm[ip][iw]/nt2; /* transpose */
#else
		mm[ip] = sf_crmul(cm[ip][iw],1.0/nt2); /* transpose */
#endif
	    }
	}
	    
	radon_set (w, xx);
	
	radon_lop (adj,false,np,nx,mm,dd); /* apply Radon */
        
	if (adj) {
	    for (ip=0; ip < np; ip++) { /* loop over slopes */
		cm[ip][iw] = mm[ip]; /* transpose */
		
	    }
	} else {
	    for (ix=0; ix < nx; ix++) { /* loop over offsets */
		cd[ix][iw] = dd[ix]; /* transpose */
	    }
	}
    } /* loop over frequencies */
    
    if (adj) {
	for (ip=0; ip < np; ip++) { /* loop over slopes */
	    /* FFT to time */
	    kiss_fftri(invs,(const kiss_fft_cpx *) cm[ip], tt);
	    for (it=0; it < nt; it++) {
		x[ip*nt+it] = tt[it];
	    }
	}
    } else { /* modeling */
	for (ix=0; ix < nx; ix++) { /* loop over offsets */
	    /* FFT to time */
	    kiss_fftri(invs,(const kiss_fft_cpx *) cd[ix], tt);
	    for (it=0; it < nt; it++) {
		y[ix*nt+it] = tt[it];
	    }
	}
    }
}

void radonoper_close(void)
/*< free allocated storage >*/
{
    free (dd);
    free (mm);
    free (*cm); free (cm);
    free (*cd); free (cd);
    free (tt);
    free (xx);
}

/* 	$Id$	 */
