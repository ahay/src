/* Computing reference slownesses */
/*
  Copyright (C) 2006 Colorado School of Mines
  Copyright (C) 2004 University of Texas at Austin
   
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

#include "slowref.h"

#include "weutil.h"
/*^*/

/*------------------------------------------------------------*/
slo3d slow_init(fslice slow_,   /* slowness slice */
		sf_axa  alx_    /* i-line (slowness/image) */,
		sf_axa  aly_    /* x-line (slowness/image) */,
		sf_axa  amz_    /* depth */,
		int    nrmax,   /* maximum number of references */
		float  dsmax,
		int   ompnth
    )
/*< initialize slowness >*/
{
    int imz, jj;

    /*------------------------------------------------------------*/
    slo3d slo;
    slo = (slo3d) sf_alloc(1,sizeof(*slo));

    slo->slow=slow_;
    slo->alx=alx_;
    slo->aly=aly_;
    slo->amz=amz_;
    slo->nrmax=nrmax;
    slo->dsmax=dsmax;
    slo->ompnth=ompnth;

    slo->ss = sf_floatalloc3(slo->alx.n,slo->aly.n,slo->ompnth);  /* slowness */
    slo->so = sf_floatalloc3(slo->alx.n,slo->aly.n,slo->ompnth);  /* slowness */
    slo->sm = sf_floatalloc2(slo->nrmax,slo->amz.n);  /* ref slowness squared */
    slo->nr = sf_intalloc              (slo->amz.n);  /* nr of ref slownesses */
    
    for (imz=0; imz<slo->amz.n; imz++) {
	fslice_get(slo->slow,imz,slo->ss[0][0]);
	
	slo->nr[imz] = slow(slo->nrmax,
			    slo->dsmax,
			    slo->alx.n*slo->aly.n,
			    slo->ss[0][0],
			    slo->sm[imz]);
    }
    for (imz=0; imz<slo->amz.n-1; imz++) {
	for (jj=0; jj<slo->nr[imz]; jj++) {
	    slo->sm[imz][jj] = 0.5*(slo->sm[imz][jj]+slo->sm[imz+1][jj]);
	}
    }

    return slo;
}

/*------------------------------------------------------------*/
int slow(int nr           /* maximum number of references */, 
	 float ds         /* minimum slowness separation */, 
	 int ns           /* number of slownesses */, 
	 const float* ss  /* [ns] slowness array */, 
	 float* sr        /* [nr] reference slownesses squared */) 
/*< compute reference slownesses, return their number >*/
{
    int is,jr,ir;
    float smin, smax, s, s2=0., qr, *ss2;

    ss2 = sf_floatalloc(ns);
    for (is=0; is<ns; is++) {
	ss2[is]=ss[is];
    }
    
    smax = sf_quantile(ns-1,ns,ss2);
    smin = sf_quantile(   0,ns,ss2);
    nr = SF_MIN(nr,1+(smax-smin)/ds);
    
    jr=0;
    for (ir=0; ir<nr; ir++) {
	qr = (ir+1.0)/nr - 0.5 * 1./nr;
	s = sf_quantile(qr*ns,ns,ss2);
	if (0==ir || SF_ABS(s-s2) > ds) {
	    sr [jr] = s*s;
	    s2 = s;
	    jr++;
	}
    }
    
    free(ss2);
    return jr;
}

/*------------------------------------------------------------*/
void slow_close( slo3d slo)
/*< close slowness >*/
{
    free(**slo->ss); free( *slo->ss); free( slo->ss);
    free(**slo->so); free( *slo->so); free( slo->so);
    ;                free( *slo->sm); free( slo->sm);
    ;                                 free( slo->nr);
}
