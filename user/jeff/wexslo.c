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

#include "wexutl.h"
/*^*/

#include "wexslo.h"

#define SOOP(a) for(ily=0;ily<cub->aly.n;ily++){ \
                for(ilx=0;ilx<cub->alx.n;ilx++){ \
		    {a} \
		}}

/*------------------------------------------------------------*/
wexslo3d wexslo_init(wexcub3d   cub,
		     sf_file  file_,   /* slowness slice */
		     int     nrmax,   /* maximum number of references */
		     float   dsmax
    )
/*< initialize slowness >*/
{
    int iz, jj;

    /*------------------------------------------------------------*/
    wexslo3d slo;
    slo = (wexslo3d) sf_alloc(1,sizeof(*slo));

    slo->file=file_;
    slo->nrmax=nrmax;
    slo->dsmax=dsmax;

    slo->s  = sf_floatalloc3(cub->alx.n,cub->aly.n,cub->az.n);
    slo->sm = sf_floatalloc2(slo->nrmax,cub->az.n);  /* ref slowness squared */
    slo->nr = sf_intalloc              (cub->az.n);  /* nr of ref slownesses */
    
    sf_floatread(slo->s[0][0],cub->alx.n*cub->aly.n*cub->az.n,slo->file);

    for (iz=0; iz<cub->az.n; iz++) {
	slo->nr[iz] = wexslo(slo->nrmax,
			      slo->dsmax,
			      cub->alx.n*cub->aly.n,
			      slo->s[iz][0],
			      slo->sm[iz]);
    }
    for (iz=0; iz<cub->az.n-1; iz++) {
	for (jj=0; jj<slo->nr[iz]; jj++) {
	    slo->sm[iz][jj] = 0.5*(slo->sm[iz][jj]+slo->sm[iz+1][jj]);
	}
    }
    
    return slo;
}

/*------------------------------------------------------------*/
int wexslo(int nr           /* maximum number of references */, 
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
void wexslo_close(wexslo3d slo)
/*< close slowness >*/
{
    free(**slo->s);  free( *slo->s);  free( slo->s);
    ;                free( *slo->sm); free( slo->sm);
    ;                                 free( slo->nr);
}

