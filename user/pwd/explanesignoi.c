/* Signal-noise separation with frequency and dip, 2-D */
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

#include "explanesignoi.h"
#include "expont.h"
#include "allp2.h"

static int n1,n2,n12;
static float eps, *a, *b, *c, *d, *tmp, *tmp2;
static allpas2 noi, sig;

void explanesignoi_init (int m1,int m2 /* data size */, 
			 float eps1    /* signal/noise scaling */, 
			 float **aa    /* frequency filter [4][m1*m2] */, 
			 int nw        /* dip filter order */, 
			 int nj1       /* dip filter step for noise */, 
			 int nj2       /* dip filter step for signal */, 
			 bool drift    /* if shift filter */,
			 float **nn    /* noise dip [m1*m2] */, 
			 float **ss    /* signal dip [m1*m2] */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    n12 = n1*n2;
    eps = eps1;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];

    noi = allpass2_init(nw,nj1,n1,n2,drift,nn);
    sig = allpass2_init(nw,nj2,n1,n2,drift,nn);

    tmp = sf_floatalloc(n12);
    tmp2 = sf_floatalloc(n12);
}

void explanesignoi_close(void)
/*< free allocated storage >*/
{
    free (tmp);
    free (tmp2);
}

void explanesignoi_lop (bool adj, bool add, int ns, int nd, 
			float *ss, float *dat)
/*< linear operator for inversion >*/
{
    int i;

    if (2*ns != nd) sf_error("%s: wrong size: 2*%d != %d",__FILE__,ns,nd);

    sf_adjnull(adj,add,ns,nd,ss,dat);

    allpass22_init(noi);
    expont_init (n1, n2, a, b);
    sf_chain(allpass21_lop, expont_lop, adj, true, ns, ns, ns, ss, dat, tmp);

    allpass22_init(sig);
    expont_init (n1, n2, c, d);

    for (i=0; i < n12; i++) {
	tmp2[i] = adj? dat[i+n12] * eps: ss[i] * eps;
    }

    if (adj) {
	sf_chain(allpass21_lop, expont_lop, 
		 true, true, ns, ns, ns, ss, tmp2, tmp);
    } else {
	sf_chain(allpass21_lop, expont_lop,
		 false, true, ns, ns, ns, tmp2, dat+n12, tmp);
    }
}
