/* 3-D SSR MVA */

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
/*^*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include "zomva3.h"
#include "zomig3.h"
#include "taper3.h"
#include "slow3.h"
#include "ssr3.h"
#include "lsr3.h"

#include "slice.h"
/*^*/

#include "weutil.h"
/*^*/
#define LOOP(a) for(imy=0;imy<amy.n;imy++){ \
	        for(imx=0;imx<amx.n;imx++){ \
		    {a} }} /* loop in x-domain */

static sf_axa amz,aw,alx,aly,amx,amy;
static bool verb;
static float eps;
static float twoway;

static float         **sm; /* reference slowness squared */
static float         **ss; /* slowness */
static float         **so; /* slowness */
static int            *nr; /* number of references */
static fslice       Bslow; /* slowness slice */
static fslice       Bwfld; /* wavefield slice */

static sf_complex **bw; /* wavefield */

static sf_complex **dw; /* wavefield */
static sf_complex **pw;
static sf_complex **pwsum;

static sf_complex **ds; /* slowness */
static sf_complex **ps;
static sf_complex **pssum;

/*------------------------------------------------------------*/
cub3d zomva3_cube(bool    verb_,
		  sf_axis amx_,
		  sf_axis amy_,
		  sf_axis amz_,
		  sf_axis alx_,
		  sf_axis aly_,
		  sf_axis aw_,
		  sf_axis ae_,
		  float   eps_,
		  int     ompnth_,
		  int     ompchunk_
    )
/*< initialize SR migration space >*/
{
    cub3d cub;
    cub = (cub3d) sf_alloc(1,sizeof(*cub));

    cub->verb=verb_;
    cub->amx = sf_nod(amx_);
    cub->amy = sf_nod(amy_);
    cub->amz = sf_nod(amz_);

    cub->alx = sf_nod(alx_);
    cub->aly = sf_nod(aly_);
    cub->ae  = sf_nod(ae_);

    cub->aw  = sf_nod(aw_);
    cub->aw.d *= 2.*SF_PI; /* from hertz to radians */
    cub->aw.o *= 2.*SF_PI;

    cub->eps      = eps_;

    cub->ompnth   = ompnth_;
    cub->ompchunk = ompchunk_;

    return cub;
}

/*------------------------------------------------------------*/
ssroperator3d zomva3_init(cub3d cub)
/*< initialize >*/
{
    ssroperator3d weop;
    weop = (ssroperator3d) sf_alloc(1,sizeof(*weop));

    weop->ww = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    weop->qq = sf_floatalloc2  (cub->amx.n,cub->amy.n);

    return weop;
}

/*------------------------------------------------------------*/
void zomva3(ssrmvaoperator3d weop,
	    cub3d cub,
	    ssr3d ssr,
	    lsr3d lsr,
	    tap3d tap,
	    slo3d slo,
	    bool  inv    /* forward/adjoint flag */, 
	    fslice Bslow,
	    fslice Pslow /* slowness perturbation [nz][nmy][nmx] */,
	    fslice Pimag /*    image perturbation [nz][nmy][nmx] */)
/*< Apply forward/adjoint ZO MVA >*/
{
    int imz,iw,imy,imx,ilx,ily;
    sf_complex w;
    int ompith=0;
    
    if(inv) {
	LOOP( weop->ps[imy][imx][ompith] = sf_cmplx(0.0,0.0); );
	for (imz=0; imz<amz.n; imz++) {
	    fslice_put(Pslow,imz,weop->ps[0]);
	}
    } else {
	LOOP( weop->pwsum[imy][imx][ompith] = sf_cmplx(0.0,0.0); );
	for (imz=0; imz<amz.n; imz++) {
	    fslice_put(Pimag,imz,weop->pwsum[0]);
	}
    }

    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if (verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);
	
	LOOP( weop->dw[imy][imx][ompith]=sf_cmplx(0.,0.); );

	if (inv) { /* adjoint: image -> slowness */
	    w = sf_cmplx(cub->eps*cub->aw.d,cub->aw.o+iw*cub->aw.d); /* causal */
	    
	    for (imz=cub->amz.n-1; imz>0; imz--) {

		/* background */
		fslice_get(slo->slice,imz,slo->so[ompith][0]);	    
		slow3_twoway(cub,slo,slo->so,ompith); /* 2way time */

		fslice_get(Bwfld,iw*amz.n+imz,weop->bw[ompith][0]);
		
		if(imz>0) { /* continuation */
		    fslice_get(slo->slice,imz-1,slo->ss[ompith][0]);
		    slow3_twoway(cub,slo,slo->ss,ompith); /* 2way time */

		    ssr3_ssf(w,weop->dw[ompith],cub,ssr,tap,slo,imz,ompith);
		    slow3_advance(cub,slo,ompith);

		} /* end continuation */
		
		/* scattering dI -> dS */
		fslice_get(Pimag,imz,weop->pwsum[ompith][0]);
		
#ifdef SF_HAS_COMPLEX_H
		LOOP(weop->dw[ompith][imy][imx] += weop->pwsum[ompith][imy][imx]; );
#else
		LOOP(weop->dw[ompith][imy][imx] = sf_cadd(weop->dw   [ompith][imy][imx],
							  weop->pwsum[ompith][imy][imx]); );
#endif
		lsr3_w2s(cub,lsr,w,bw,so,dw,ps);
		
		fslice_get(Pslow,imz,weop->pssum[ompith][0]);
#ifdef SF_HAS_COMPLEX_H
		LOOP(weop->pssum[ompith][imy][imx] += weop->ps[ompith][imy][imx];);
#else
		LOOP(weop->pssum[ompith][imy][imx] = sf_cadd(weop->pssum[ompith][imy][imx],
							     weop->ps   [ompith][imy][imx]););
#endif
		fslice_put(Pslow,imz,weop->pssum[ompith][0]);
		/* end scattering */

	    }
	} else {
	    w = sf_cmplx(cub->eps*cub->aw.d,-(cub->aw.o+iw*cub->aw.d)); /* anti-causal */
	    
	    for (imz=0; imz<cub->amz.n-1; imz++) {
		
	    }
	}
	
}

/*------------------------------------------------------------*/
void zomva3_close(ssrmvaoperator3d weop)
/*< free allocated storage >*/
{
    free( **weop->bw);free( *weop->bw); free( weop->bw);

    free( **weop->dw);free( *weop->dw); free( weop->dw);
    free( **weop->pw);free( *weop->pw); free( weop->pw);

    free( **weop->ds);free( *weop->ds); free( weop->ds);
    free( **weop->ps);free( *weop->ps); free( weop->ps);

    free( **weop->pwsum);free( *weop->pwsum); free( weop->pwsum);
    free( **weop->pssum);free( *weop->pssum); free( weop->pssum);
}

