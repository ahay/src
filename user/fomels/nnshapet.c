/* 2-D irregular data interpolation using nearest neighbors and shaping regularization. */
/*
  Copyright (C) 2013 University of Texas at Austin

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

#include "distance.h"

static bool sym;
static int nd,nt, n1, n2, n12, n123;
static float o1,d1, o2,d2;
static float **d, **w, *m, *fold, *bin, *zero, eps;
static sf_triangle tr1, tr2, tr3;
static sf_upgrad upg;

void nnshapet_init(bool sym1, int nd1, int nt1, int n11, int n21, 
		  float o11, float o21, float d11, float d21,
		  int rect1, int rect2, int rect3, int nw, int order,
		   float **xy, float **weight, float eps)
/*< initialize >*/
{
    int i, id, n[2], *pp;
    float dd[2], **pts;
    
    sym = sym1;
    nd = nd1;
    nt = nt1; 

    n1 = n11;
    n2 = n21;

    n12  = n1*n2;
    n123 = n12*nt;

    o1 = o11; d1=d11;
    o2 = o21; d2=d21;

    n[0]=n1; dd[0]=d1;
    n[1]=n2; dd[1]=d2;

    tr1 = (rect1 > 1)? sf_triangle_init (rect1,n1,false): NULL;
    tr2 = (rect2 > 1)? sf_triangle_init (rect2,n2,false): NULL;
    tr3 = (rect3 > 1)? sf_triangle_init (rect3,nt,false): NULL;

    sf_int2_init (xy, o1, o2, d1, d2, n1, n2, sf_lg_int, nw, nd);

    d = sf_floatalloc2(nd,nt);
    m = sf_floatalloc(n123);
    w = weight;

    fold = sf_floatalloc(n12);
    bin = sf_floatalloc(n12);
    zero = sf_floatalloc(n12);
    
    /* compute the fold for normalization */
    for (id=0; id<nd; id++) {
	d[0][id]=1.0f;
    }
    sf_int2_lop (true,false,n12,nd,fold,d[0]);

    upg = sf_upgrad_init(2,n,dd);

    /* 1. find distance */
    for(i = 0; i < n12; i++) {
	bin[i] = 1.0f;
	zero[i] = 0.0f;
    }
    
    distance_init (1,n2,n1,nd);

    pts = sf_floatalloc2 (3,nd);
    pp = sf_intalloc (n12);
    
    for (id=0; id<nd; id++) {
	pts[id][0] = xy[id][0];
	pts[id][1] = xy[id][1];
	pts[id][2] = 0.0f;
    }
    
    distance(nd,pts,m,bin,pp,
	     1,n2,n1,
	     0.0f,o2,o1,
	     1.0f,d2,d1,
	     order);
    sf_upgrad_set(upg,m);

    free(*pts);
    free(pts);
    free(pp);
}

void nnshapet_close(void)
/*< free allocated storage >*/
{
    if (NULL != tr1) sf_triangle_close(tr1);
    if (NULL != tr2) sf_triangle_close(tr2);
    if (NULL != tr3) sf_triangle_close(tr3);

    sf_int2_close();
    free(d);
    free(m);
    free(fold);
    free(bin);
    free(zero);
}



void nnshapet_smooth(float *modl)
/*< shaping operator (triangle smoothing) >*/
{
    int i1, i2, it;
    
    if (NULL != tr1) {
	for (it=0; it < nt; it++) {	
	   for (i2=0; i2 < n2; i2++) {
	       sf_smooth2 (tr1, 0, 1, false, modl+n1*(i2+it*n2));
	   }
        }
    }
    if (NULL != tr2) {
        for (it=0; it < nt; it++) {	
	    for (i1=0; i1 < n1; i1++) {
	        sf_smooth2 (tr2, i1, n1, false, modl+it*n1*n2);
            }
	}
    }   
    if (NULL != tr2) {
        for (i2=0; i2 < n2; i2++) {	
	    for (i1=0; i1 < n1; i1++) {
	        sf_smooth2 (tr2, i1 + i2*n1,n1*n2, false, modl);
            }
	}
    } 
}    

void nnshapet_back(float **data, float *modl)
/*< backward operator (voronoi diagrams) >*/
{
    int i, it, id;
    
    for (it=0; it < nt; it++) {
	if (NULL != w) {
	    for (id=0; id < nd; id++) {
		data[it][id] *= w[it][id]/(w[it][id]*w[it][id]+eps);
	    }
	}
	sf_int2_lop (true,false,n12,nd,bin,data[it]);
	for (i=0; i < n12; i++) {
	    /* normalize by the fold */
	    if (fold[i] > FLT_EPSILON) bin[i] /=fold[i];	
	}
	sf_upgrad_solve(upg,zero,modl+it*n12,bin);
    } 
}

void nnshapet_forw(const float *modl, float **data)
/*< forward operator (grid interpolation) >*/
{
    int it, id;

    for (it=0; it < nt; it++) {
       sf_int2_lop (false,false,n12,nd,(float*) (modl+it*n12),data[it]);
       if (NULL != w) {
	   for (id=0; id < nd; id++) {
	       data[it][id] *= w[it][id];
	   }
       }
    }
}

void nnshapet(int n, const float *inp, float *out, void *data)
/*< I + S (BF - I) >*/
{
    int i;

    if (sym) {
	/* (B F - I) S */
	for (i=0; i < n; i++) {
	    m[i] = inp[i];
	}
	nnshapet_smooth(m);

	nnshapet_forw(m,d);
	nnshapet_back(d,out);

	for (i=0; i < n; i++) {
	    out[i] -= m[i];
	}
    } else {
	/* B F - I  */

	nnshapet_forw(inp,d);
	nnshapet_back(d,out);

	for (i=0; i < n; i++) {
	    out[i] -= inp[i];
	}
    }

    /* S */
    nnshapet_smooth(out);

    /* + I */
    for (i=0; i < n; i++) {
	out[i] += inp[i];
    }
}
