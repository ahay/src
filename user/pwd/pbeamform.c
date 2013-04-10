/* Beam forming with plane-wave construction. */
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

#include "pbeamform.h"
#include "pwdsl.h"

static bool gauss;
static int n1, n2, rect;
static float *t, *t2;

void pbeamform_init(bool gauss1    /* pseudo-gaussian */,
		    int m1, int m2 /* data dimensions */, 
		    int order, 
		    int rect1      /* triangle radius */)
/*< initialize >*/
{
    const float eps=0.01;

    n1=m1;
    n2=m2;
    rect=rect1;
    gauss=gauss1;

    t = sf_floatalloc(n1*n2);
    if (gauss) t2 = sf_floatalloc(n1*n2);

    pwdsl_init(n1,n2,order,1,rect,eps);
}

void pbeamform_set(float** p /* slope [m2][m1] */)
/*< set slopes >*/
{
    pwdsl_set(p);
}

void pbeamform_close(void)
/*< free allocated storage >*/
{
    free(t);
    if (gauss) free(t2);
    pwdsl_close();
}

void pbeamform_lop(bool adj, bool add, int nc, int nd, float* c, float* d)
/*< linear operator >*/
{
    int i1, ic, id;

    if (nd != n1*n2) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj,add,nc,nd,c,d);

    if (adj) {
	if (gauss) {
	    pwdsl_lop(false,false,nd,nd,d,t2);
	    pwdsl_lop(true,false,nd,nd,t,t2);
	} else {
	    pwdsl_lop(true,false,nd,nd,t,d);
	}
	for (ic=id=0; id < n2; id++) {
	    if (0==id%rect) {
		for (i1=0; i1 < n1; i1++) {
		    c[ic*n1+i1] += t[id*n1+i1];
		}
		ic++;
	    } 
	}
    } else {
	for (ic=id=0; id < n2; id++) {
	    if (0==id%rect) {
		for (i1=0; i1 < n1; i1++) {
		    t[id*n1+i1] = c[ic*n1+i1];
		}
		ic++;
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    t[id*n1+i1] = 0.;
		}
	    }
	}
	if (gauss) {
	    pwdsl_lop(false,false,nd,nd,t,t2);
	    pwdsl_lop(true,true,nd,nd,d,t2);
	} else {
	    pwdsl_lop(false,true,nd,nd,t,d);
	}
    }
}


