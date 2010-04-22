/* First-arrival tomography. */
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

#include "upgrad.h"
#include "fatomo.h"

static int ns, nm;
static upgrad *upg;
static float *t;
static bool **mask;

void fatomo_init(int ns1             /* number of sources */,
		 int dim             /* model dimensions */,
		 const int *n        /* model size */,
		 const float *d      /* model sampling */,
		 bool **masks  /* [ns1][nm1] receiver mask */)
/*< initialize >*/
{
    int is, i;

    ns = ns1;

    upg = (upgrad*) sf_alloc(ns,sizeof(upgrad));

    for (is=0; is < ns; is++) {
	upg[is] = upgrad_init(dim,n,d);
    }

    nm = 1;
    for (i=0; i < dim; i++) {
	nm *= n[i];
    }

    t = sf_floatalloc(nm);

    mask = masks;
}

void fatomo_close(void)
/*< free allocated storage >*/
{
    int is;

    for (is=0; is < ns; is++) {
	upgrad_close(upg[is]);
    }
    free(upg);
    free(t);
}

void fatomo_set(const float **times /* [ns1][nm] background times */)
/*< set background traveltime >*/
{ 
    int is;

    for (is=0; is < ns; is++) {
	upgrad_set(upg[is],times[is]);
    }
}

void fatomo_lop(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int ir, is, im;

    if (nx != nm) sf_error("%s: Wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,nr,x,r);

    ir = 0;
    for (is=0; is < ns; is++) {
	if (adj) {
	    for (im=0; im < nm; im++) {
		if (mask[is][im]) {
		    t[im] = r[ir];
		    ir++;
		} else {
		    t[im] = 0.;
		}
	    }
	    upgrad_inverse(upg[is],true,x,t,NULL);
	} else {
	    upgrad_solve(upg[is],x,t,NULL);
	    for (im=0; im < nm; im++) {
		if (mask[is][im]) {
		    r[ir] += t[im];
		    ir++;
		}
	    }
	}
    }
}
    

