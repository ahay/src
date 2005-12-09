/* one-way Rienmannian Wavefield Extrapolation 
   zero-offset migration subroutines
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

#include "rwezom.h"

static axa ag,at,aw,ar;
static int method;
static bool verb;
static char *met;

void rwezom_init(
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

void rwezom_main(
    bool            adj,
    complex float **wfl,
    float         **img,
    float         ** aa,
    float         ** bb,
    float         ** mm,
    float         ** a0,
    float         ** b0)
/*< run one-way extrapolation >*/
{
    int iw,it;
    float w;

    switch(method) {
	case 3: rweone_psc_coef(aa,bb,a0,b0); break;
	case 2: rweone_ffd_coef(aa,bb,a0,b0); break;
	case 1: ; /* SSF */                   break;
	case 0: rweone_xfd_coef(aa,bb);       break;
    }
    
    for(iw=0;iw<aw.n;iw++) {
	w=aw.o+iw*aw.d;
	sf_warning("%s %d %d",met,iw,aw.n);

	if(adj) {  /* modeling */
	    w*=-2; /*      causal, two-way time */
	    
	    it=at.n-1; rweone_zoi(adj,wfl[iw],img[it]);
	    for(it=at.n-2;it>=0;it--) {
		if(method!=0) rweone_fk(w,wfl[iw],aa[it],a0[it],b0[it],mm[it],it);
		else          rweone_fx(w,wfl[iw],aa[it],it);
		rweone_zoi(adj,wfl[iw],img[it]);
	    }
	} else {   /* migration */
	    w*=+2; /* anti-causal, two-way time */
	    
	    for(it=0;it<=at.n-2;it++) {
		rweone_zoi(adj,wfl[iw],img[it]);		
		if(method!=0) rweone_fk(w,wfl[iw],aa[it],a0[it],b0[it],mm[it],it);
		else          rweone_fx(w,wfl[iw],aa[it],it);
	    }
	    it=at.n-1; rweone_zoi(adj,wfl[iw],img[it]);
	}
    } /* w loop */
}

/*------------------------------------------------------------*/

