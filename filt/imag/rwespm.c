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

#include "rweone.h"
/*^*/

#include "rwespm.h"

static axa ag,at,aw,ar;
static int method;
static bool verb;
static char *met;

void rwespm_init(
    axa ag_,
    axa at_,
    axa aw_,
    axa ar_,
    int method_,
    bool verb_)
/*< initialize >*/
{
    ag = ag_;
    at = at_;
    aw = aw_;
    ar = ar_;

    method = method_;
    verb   = verb_;

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    rweone_init(ag_,at_,aw_,ar_,method_,verb_);

    switch (method) {
	case 3: met="PSC"; break;
	case 2: met="FFD"; break;
	case 1: met="SSF"; break;
	case 0: met="XFD"; break;
    }
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
	case 3: rweone_psc_coef(aa,bb,a0,b0); break;
	case 2: rweone_ffd_coef(aa,bb,a0,b0); break;
	case 1: ;/* SSF */                    break;
	case 0: rweone_xfd_coef(aa,bb);       break;
    }

    for(iw=0;iw<aw.n;iw++) {
	w=aw.o+iw*aw.d;
	sf_warning("%s %d %d",met,iw,aw.n);

	ws = -w; /*      causal */
	wr = +w; /* anti-causal */
	for(it=0;it<=at.n-2;it++) {
	    rweone_spi(swf[iw],rwf[iw],img[it]);

	    if(method!=0) {
		rweone_fk(ws,swf[iw],aa[it],a0[it],b0[it],mm[it],it);
		rweone_fk(wr,rwf[iw],aa[it],a0[it],b0[it],mm[it],it);
	    } else {
		rweone_fx(ws,swf[iw],aa[it],it);
		rweone_fx(wr,rwf[iw],aa[it],it);
	    }
	}
	it=at.n-1; rweone_spi(swf[iw],rwf[iw],img[it]);
    }
}

/*------------------------------------------------------------*/

