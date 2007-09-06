/* Computing reference slownesses */
/*
  Copyright (C) 2007 Colorado School of Mines
   
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

#include "slow3.h"

#include "weutil.h"
/*^*/

#define SOOP(a) for(ily=0;ily<cub->aly.n;ily++){ \
                for(ilx=0;ilx<cub->alx.n;ilx++){ \
		    {a} \
		}}

/*------------------------------------------------------------*/
slo3d slow3_init(cub3d   cub,
		 fslice  slice_,   /* slowness slice */
		 int     nrmax,   /* maximum number of references */
		 float   dsmax,
		 float  twoway
    )
/*< initialize slowness >*/
{
    int imz, jj;
    int ompith=0;

    /*------------------------------------------------------------*/
    slo3d slo;
    slo = (slo3d) sf_alloc(1,sizeof(*slo));

    slo->slice=slice_;
    slo->nrmax=nrmax;
    slo->dsmax=dsmax;

    if(twoway) { slo->twoway = 2;
    } else {     slo->twoway = 1;
    }

    slo->ss = sf_floatalloc3(cub->alx.n,cub->aly.n,cub->ompnth);  /* slowness */
    slo->so = sf_floatalloc3(cub->alx.n,cub->aly.n,cub->ompnth);  /* slowness */
    slo->sm = sf_floatalloc2(slo->nrmax,cub->amz.n);  /* ref slowness squared */
    slo->nr = sf_intalloc              (cub->amz.n);  /* nr of ref slownesses */
    
    for (imz=0; imz<cub->amz.n; imz++) {
	fslice_get(slo->slice,imz,slo->ss[0][0]);
	slow3_twoway(cub,slo,slo->ss,ompith);

	slo->nr[imz] = slow3(slo->nrmax,
			     slo->dsmax,
			     cub->alx.n*cub->aly.n,
			     slo->ss[0][0],
			     slo->sm[imz]);
    }
    for (imz=0; imz<cub->amz.n-1; imz++) {
	for (jj=0; jj<slo->nr[imz]; jj++) {
	    slo->sm[imz][jj] = 0.5*(slo->sm[imz][jj]+slo->sm[imz+1][jj]);
	}
    }

    return slo;
}

/*------------------------------------------------------------*/
int slow3(int nr           /* maximum number of references */, 
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
void slow3_close( slo3d slo)
/*< close slowness >*/
{
    free(**slo->ss); free( *slo->ss); free( slo->ss);
    free(**slo->so); free( *slo->so); free( slo->so);
    ;                free( *slo->sm); free( slo->sm);
    ;                                 free( slo->nr);
}

/*------------------------------------------------------------*/
void slow3_advance( cub3d cub,
		    slo3d slo,
		    int ompith)
/*< close slowness >*/
{
    int ilx,ily;

    SOOP( slo->so[ompith][ily][ilx] = slo->ss[ompith][ily][ilx]; );

}

/*------------------------------------------------------------*/
void slow3_twoway( cub3d     cub,
		   slo3d     slo,
		   float ***slow,
		   int    ompith)
/*< close slowness >*/
{
    int ilx,ily;

    SOOP( slow[ompith][ily][ilx]*= slo->twoway; );

}


