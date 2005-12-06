/* one-way Rienmannian Wavefield Extrapolation 
   shot-record migration subroutines
*/
/*
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

#include <math.h>

#include <rsf.h>
/*^*/

#include "rwespm.h"

#include "rweone.h"
/*^*/

static axa ag,at,aw,ar;
static int method;

void rwespm_init(
    axa ag_,
    axa at_,
    axa aw_,
    axa ar_,
    int method_)
/*< initialize >*/
{
    ag = ag_;
    at = at_;
    aw = aw_;
    ar = ar_;

    method =method_;
    
    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    rweone_init(ag_,at_,aw_,ar_,method_);
}

void rwespm_main(
    complex float **swf,
    complex float **rwf,
    float         **img,
    float         ** aa,
    float         ** bb,
    float         ** mm,
    float         ** a0,
    float         ** b0)
/*< run one-way extrapolation >*/
{
    int iw,it;
    float w, ws,wr;

    switch(method) {
	case 3: /* PSC */
	    sf_warning("compute PSC coef");
	    rweone_psc_coef(aa,bb,a0,b0);
	    break;
	case 2: /* FFD */
	    sf_warning("compute FFD coef");
	    rweone_ffd_coef(aa,bb,a0,b0);
	    break;
	case 1: /* SSF */
	    break;
	case 0: /* XFD */
	    sf_warning("compute XFD coef");
	    rweone_xfd_coef(aa,bb);
	    break;
    }

    switch(method) {
	case 3:
	case 2:
	case 1:
	    for(iw=0;iw<aw.n;iw++) {
		w=aw.o+iw*aw.d;

		if     (method==3) { sf_warning("PSC %d %d",iw,aw.n); }
		else if(method==2) { sf_warning("FFD %d %d",iw,aw.n); }
		else               { sf_warning("SSF %d %d",iw,aw.n); }

		ws = -w; /*      causal */
		wr = +w; /* anti-causal */
		for(it=0;it<=at.n-2;it++) {
		    rwespm_img(swf[iw],rwf[iw],img[it]);

		    rweone_fk(ws,swf[iw],aa[it],a0[it],b0[it],mm[it],it);
		    rweone_fk(wr,rwf[iw],aa[it],a0[it],b0[it],mm[it],it);
		}
		it=at.n-1; rwespm_img(swf[iw],rwf[iw],img[it]);
	    }
	    break;
	case 0:
	default:
	    for(iw=0;iw<aw.n;iw++) {
		w=aw.o+iw*aw.d; 

		sf_warning("XFD %d %d",iw,aw.n);

		ws = -w; /*      causal */
		wr = +w; /* anti-causal */
		for(it=0;it<=at.n-2;it++) {
		    rwespm_img(swf[iw],rwf[iw],img[it]);

		    rweone_fx(ws,swf[iw],aa[it],it);
		    rweone_fx(wr,rwf[iw],aa[it],it);
		}
		it=at.n-1; rwespm_img(swf[iw],rwf[iw],img[it]);
	    }
	    break;
    } /* end switch method */
}

/*------------------------------------------------------------*/

void rwespm_img(
    complex float *swf,
    complex float *rwf,
    float         *iii)
/*< imaging condition >*/
{
    int ig;

    rweone_tap(swf);
    rweone_tap(rwf);
    for(ig=0;ig<ag.n;ig++) {
	iii[ig] += crealf( conjf(swf[ig]) * rwf[ig] );
    }
}

