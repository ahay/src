/* Step of L1-like iteration. */
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

#include <rsf.h>

#include "l1.h"

static float* S;  /* model step */
static float* Ss; /* residual step */
static bool Allocated = false; /* if S and Ss are allocated */

void l1step( bool forget     /* restart flag */, 
	     int nx, int ny  /* model size, data size */, 
	     float* x        /* current model [nx] */, 
	     const float* g  /* gradient [nx] */, 
	     float* rr       /* data residual [ny] */, 
	     const float* gg /* conjugate gradient [ny] */) 
/*< Step of conjugate-gradient-like iteration. >*/
{
    float alfa, beta, l1;
    int i;

    if (!Allocated) {
	Allocated = forget = true;
	S  = sf_floatalloc (nx);
	Ss = sf_floatalloc (ny);
    }
    if (forget) {
	for (i = 0; i < nx; i++) S[i] = 0.;
	for (i = 0; i < ny; i++) Ss[i] = 0.;
	beta = 0.0;
	/* Minimize |R + alfa*G| */
	unil1(rr,gg,&alfa);
	alfa = - alfa;
    } else {
	/* Minimize |R + alfa*G + beta*S| */
	bil1(rr,gg,Ss,&alfa,&beta);
	alfa = -alfa;
	beta = -beta;
    }
    cblas_sscal(nx,beta,S,1);
    cblas_saxpy(nx,alfa,g,1,S,1);

    cblas_sscal(ny,beta,Ss,1);
    cblas_saxpy(ny,alfa,gg,1,Ss,1);

    for (i = 0; i < nx; i++) {
	x[i] +=  S[i];
    }

    l1 = 0.0;
    for (i = 0; i < ny; i++) {
	rr[i] += Ss[i];
	l1 += fabsf(rr[i]);
    }
    sf_warning("l1 step -> %g",l1);
}

void l1step_close (void) 
/*< Free allocated space. >*/ 
{
    if (Allocated) {
	free (S);
	free (Ss);
	Allocated = false;
    }
}

/* 	$Id: cgstep.c 6511 2010-08-18 23:19:11Z sfomel $	 */

