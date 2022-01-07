/* Signal and noise separation with plane-wave destruction */
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
/*^*/

#include "planesignoi.h"

#include "allp2.h"

static float eps, *tmp;
static allpas2 noi, sig;

void planesignoi_init (int nw           /* filter size */, 
		       int nj1, int nj2 /* dealiasing stretch */, 
		       int nx, int ny   /* data size */, 
		       bool drift       /* if shift filter */,
		       float **nn       /* noise slope */, 
		       float **ss       /* signal slope */, 
		       float eps1       /* regularization parameter */)
/*< initialize >*/
{
    eps = eps1;
    noi = allpass2_init(nw,nj1,nx,ny,drift,nn);
    sig = allpass2_init(nw,nj2,nx,ny,drift,ss);
    tmp = sf_floatalloc(nx*ny);
}

void planesignoi_lop (bool adj, bool add, int ns, int nd, float *s, float *d)
/*< linear operator >*/
{
    int is;

    if (nd != 2*ns) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,ns,nd,s,d);

    allpass22_init(noi);
    allpass21_lop (adj, true, ns, ns, s, d);

    allpass22_init (sig);
    if (adj) {
	for (is=0; is < ns; is++) {
	    tmp[is] = eps*d[ns+is];
	}
	allpass21_lop (true, true, ns, ns, s, tmp);
    } else {
	allpass21_lop (false, false, ns, ns, s, tmp);
	for (is=0; is < ns; is++) {
	    d[ns+is] += eps*tmp[is];
	}
    }
}

void planesignoi_close(void)
/*< free allocated storage >*/
{
    free(noi);
    free(sig);
    free(tmp);
}
