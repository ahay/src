/* Compute linear matching-radon operator */
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
#include <rsf.h>
/*^*/
#include "matchoper.h"
#include "radonoper.h"
#include "slant.h"

static float *d, pclip;
static int nt, nx, np, shift, niter;
static bool freq;

void matchoper_init (int shift1, int mdim1, int niter1,
		     int *m1, int *rect1, bool verb1,
		     int nt1, int nx1, int np1,
		     float dt1, float t01,
		     float dx1, float ox1, float x01,
		     float dp1, float p01,
		     bool par1, bool freq1, bool rho1,
		     float anti1, float p11, float pclip1) 
/*< initialize >*/
{
    nt = nt1;
    nx = nx1;
    np = np1;
    freq = freq1;
    shift = shift1;
    niter = niter1;
    pclip = pclip1;

    d = sf_floatalloc(nt1*nx1*shift1);

    sf_multidivn_init(shift1, mdim1, nt1*nx1, m1, rect1, d, NULL, verb1); 

    if (freq) {
	radonoper_init (nt1,dt1,t01,nx1,dx1,ox1,x01,np1,dp1,p01,par1);
    } else {
	slant_init (true,rho1,x01,dx1,nx1,p01,dp1,np1,t01, dt1,nt1,p11,anti1);
    }
}

void matchoper_lop (float* filter, float* dd, float* wei) 
/*< linear matching operator >*/
{
    int i, j, k, n12, nd, index, n1=0, nm;
    float mean, *g, *tmp, *atmp=NULL, t=0., td=0.;

    nd = nt*nx;
    n12 = nd*shift;
    nm = nt*np;

    g = sf_floatalloc(nd);
    tmp = sf_floatalloc(nm);
    if (NULL != wei) {
	atmp = sf_floatalloc(nm);
	n1 = 0.5+nm*(1.-0.01*pclip);
	if (n1 < 0) n1=0;
	if (n1 >= nm) n1=nm-1;
    }

    if (freq) {
	radonoper_lop (true, false, nm, nd, tmp, dd);
	if (NULL != wei) {
	    for (i=0; i < nm; i++) {
		atmp[i] = fabsf(tmp[i]);
	    }
	    t = sf_quantile(n1,nm,atmp);
	    for (i=0; i < nm; i++) {
		atmp[i] = tmp[i];
	    }
	    for (i=0; i < nm; i++) {
		td = atmp[i];
		if (td < -t) {
		    atmp[i] = td+t;
		} else if (td > t) {
		    atmp[i] = td-t;
		} else {
		    atmp[i] = 0.;
		}    
	    }
	    for (i=0; i < nm; i++) {
		wei[i] = atmp[i]/(tmp[i]+FLT_EPSILON);
	    }
	    for (i=0; i < nm; i++) {
		tmp[i] *= wei[i];
	    }
	}
	radonoper_lop (false,false, nm, nd, tmp, g);
    } else {
	slant_lop(true,false,nm,nd,tmp,dd);	
	if (NULL != wei) {
	    for (i=0; i < nm; i++) {
		atmp[i] = fabsf(tmp[i]);
	    }
	    t = sf_quantile(n1,nm,atmp);
	    for (i=0; i < nm; i++) {
		atmp[i] = tmp[i];
	    }
	    for (i=0; i < nm; i++) {
		td = atmp[i];
		if (td < -t) {
		    atmp[i] = td+t;
		} else if (td > t) {
		    atmp[i] = td-t;
		} else {
		    atmp[i] = 0.;
		}    
	    }
	    for (i=0; i < nm; i++) {
		wei[i] = atmp[i]/(tmp[i]+FLT_EPSILON);
	    }
	    for (i=0; i < nm; i++) {
		tmp[i] *= wei[i];
	    }
	}
	slant_lop(false,false,nm,nd,tmp,g);
    }

/*    for (k=-shift/2; k < (shift+1)/2; k++) {
	for (j=0; j < nx; j++) {
	    for (i=0; i < nt; i++) {
		if ((i+k) >=0 && (i+k) < nt) {
		    d[(k+shift/2)*nt*nx+j*nt+i] = g[j*nt+i+k];
		} else {
		    d[(k+shift/2)*nt*nx+j*nt+i] = 0.;
		}
	    }	
	}	
    }
*/
    for (j=0; j < nx; j++) {
	for (i=0; i < nt; i++) {
	    d[j*nt+i] = g[j*nt+i];
	}
    }
    index = 1;
    for (k=1; k < (shift+1)/2; k++) {
	for (j=0; j < nx; j++) {
	    for (i=0; i < nt; i++) {
		if ((i+k) >=0 && (i+k) < nt) {
		    d[index*nt*nx+j*nt+i] = g[j*nt+i+k];
		} else {
		    d[index*nt*nx+j*nt+i] = 0.;
		}
		if ((i-k) >=0 && (i-k) < nt) {
		    d[(index+1)*nt*nx+j*nt+i] = g[j*nt+i-k];
		} else {
		    d[(index+1)*nt*nx+j*nt+i] = 0.;
		}
	    }	
	}
	index +=2;	    
    } /* Python shifts structure */
    
    mean = 0.;
    for(i=0; i < n12; i++) {
	mean += d[i]*d[i];
    }
    if (mean == 0.) {
	for (k=0; k < n12; k++) {
	    filter[k] = d[k];
	}
	exit(0);
    }

    mean = sqrtf (mean/n12);
    
    for(i=0; i < n12; i++) {
	d[i] /= mean;
    }
    for(i=0; i < nd; i++) {
	dd[i] /= mean;
    }
    
    sf_multidivn (dd,filter,niter);    
}

/* 	$Id$	 */
