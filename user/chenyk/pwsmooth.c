/* Structure-oriented smoothing */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

#include "pwspray.h"
#include "pwsmooth.h"

static int n1, n2, n3, ns1, ns2, n12, n13, order1, order2;
static float ***u, *w, **w1, *t;
static float **p1, **p2, *smooth1, *smooth2, *smooth3;
static float ***idip, ***xdip, ***itmp, ***itmp2, ***xtmp, ***xtmp2;
static float eps;


void pwsmooth_init(int ns      /* spray radius */,
		   int m1      /* trace length */,
		   int m2      /* number of traces */,
		   int order   /* PWD order */,
		   float eps   /* regularization */)
/*< initialize >*/
{
    int is;

    n1 = m1;
    n2 = m2;
    n12 = n1*n2;

    ns2 = pwspray_init(ns,n1,n2,order,eps);

    u = sf_floatalloc3(n1,ns2,n2);
    w = sf_floatalloc(ns2);
    w1 = sf_floatalloc2(n1,n2);

    for (is=0; is < ns2; is++) {
	w[is]=ns+1-SF_ABS(is-ns);
    }

    /* Normalization */
    t = sf_floatalloc(n12);
}

void pwsmooth_set(float **dip /* local slope */)
/*< set local slope >*/
{
    int i1;
    
    pwspray_set(dip);

    for (i1=0; i1 < n12; i1++) {
	w1[0][i1]=1.0f;
    }

    pwsmooth_lop(false,false,n12,n12,w1[0],t);

    for (i1=0; i1 < n12; i1++) {
	if (0.0f != t[i1]) {
	    w1[0][i1]=1.0/t[i1];
	} else {
	    w1[0][i1]=0.0f;
	}
    }
}

void pwsmooth_close(void)
/*< free allocated storage >*/
{
    free(**u);
    free(*u);
    free(u);
    free(w);
    free(*w1);
    free(w1);
    free(t);
    pwspray_close();
}

void pwsmooth_lop(bool adj, bool add, 
		  int nin, int nout, float* trace, float *smooth)
/*< linear operator >*/
{
    int i1, i2, is;
    float ws;

    if (nin != nout || nin != n1*n2) 
	sf_error("%s: wrong size %d != %d",__FILE__,nin,nout);
    
    sf_adjnull(adj,add,nin,nout,trace,smooth);

    if (adj) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		ws=w1[i2][i1]; 
		for (is=0; is < ns2; is++) {
		    u[i2][is][i1] = smooth[i2*n1+i1]*w[is]*ws;
		}
	    }
	}

	pwspray_lop(true,  true,  nin, nin*ns2, trace, u[0][0]);
    } else {
	pwspray_lop(false, false, nin, nin*ns2, trace, u[0][0]);

	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		ws=w1[i2][i1]; 
		for (is=0; is < ns2; is++) {
		    smooth[i2*n1+i1] += u[i2][is][i1]*w[is]*ws;
		}
	    }
	}
    }
}





void pwsmooth2_init(int ns      /* spray radius */,  
		    int m1      /* trace length */,
		    int m2      /* number of traces */,
		    int order   /* PWD order */,
		    float eps   /* regularization */)
/*< initialize >*/
{
    n1=m1;
    n2=m2;
    n12=n1*n2;

    smooth1 = sf_floatalloc(n12);
    smooth2 = sf_floatalloc(n12);
    smooth3 = sf_floatalloc(n12);
    pwsmooth_init(ns,n1,n2,order,eps);
}

void pwsmooth2_set(float **dip1 /* local slope 1 */,
		   float **dip2 /* local slope 2 */)
/*< set local slope >*/
{
    p1 = dip1;
    p2 = dip2;
}

void pwsmooth2_close(void)
/*< free allocated storage >*/
{
    free(smooth1);
    free(smooth2);
    free(smooth3);

    pwsmooth_close();
}

void pwsmooth2_lop(bool adj, bool add, 
		  int nin, int nout, float* trace, float *smooth)
/*< linear operator >*/
{
    int i;

    if (nin != nout || nin != n12) 
	sf_error("%s: wrong size %d != %d",__FILE__,nin,nout);
    
    sf_adjnull(adj,add,nin,nout,trace,smooth);

    if (adj) {
	pwsmooth_set(p1);
	pwsmooth_lop(true,false,nin,nout,smooth2,smooth);
	pwsmooth_set(p2);
	pwsmooth_lop(true,false,nin,nout,smooth1,smooth);
	pwsmooth_lop(true,false,nin,nout,smooth3,smooth2);
	for (i=0; i < n12; i++) {
	    trace[i] += smooth2[i]-smooth3[i];
	}
	pwsmooth_set(p1);
	pwsmooth_lop(true,false,nin,nout,smooth3,smooth1);
	for (i=0; i < n12; i++) {
	    trace[i] += smooth1[i]-smooth3[i];
	}
    } else {
	pwsmooth_set(p1);
	pwsmooth_lop(false,false,nin,nout,trace,smooth1);
	for (i=0; i < n12; i++) {
	    smooth1[i] = trace[i]-smooth1[i];
	}
	pwsmooth_set(p2);
	pwsmooth_lop(false,false,nin,nout,trace,smooth2);
	for (i=0; i < n12; i++) {
	    smooth2[i] = trace[i]-smooth2[i];
	}
	pwsmooth_lop(false,true,nin,nout,smooth1,smooth);
	pwsmooth_set(p1);
	pwsmooth_lop(false,true,nin,nout,smooth2,smooth);
    }
}




void pwsmooth3_init(int ns1_in      /* spray radius */,
		    int ns2_in      /* spray radius */,
		    int m1          /* trace length */,
		    int m2          /* number of traces */,
		    int m3          /* number of traces */,
		    int order1_in   /* PWD order */,
		    int order2_in   /* PWD order */,
		    float eps_in    /* regularization */,
		    float ****dip   /* local slope */)
/*< initialize >*/
{
    int i2, i3;

    n1 = m1;
    n2 = m2;
    n3 = m3;
    n12 = n1*n2;
    n13 = n1*n3;

    ns1 = ns1_in;
    ns2 = ns2_in;

    order1 = order1_in;
    order2 = order2_in;
    
    eps = eps_in;

    idip = dip[0];
    xdip = (float***) sf_alloc(n2,sizeof(float**));
    for (i2=0; i2 < n2; i2++) {
	xdip[i2] = (float**) sf_alloc(n3,sizeof(float*));
	for (i3=0; i3 < n3; i3++) {
	    xdip[i2][i3] = dip[1][i3][i2];
	}
    }
    
    itmp = sf_floatalloc3(n1,n2,n3);
    xtmp = sf_floatalloc3(n1,n3,n2);
    itmp2 = sf_floatalloc3(n1,n3,n2);

//     xtmp2 = (float***) sf_alloc(n3,sizeof(float**));
//     for (i3=0; i3 < n3; i3++) {
// 	xtmp2[i3] = (float**) sf_alloc(n2,sizeof(float*));
// 	for (i2=0; i2 < n2; i2++) {
// 	    xtmp2[i3][i2] = xtmp[i2][i3];
// 	}
//     }
}

void pwsmooth3_close(void)
/*< free allocated storage >*/
{
    int i2, i3;

    for (i2=0; i2 < n2; i2++) {
	free(xdip[i2]);
    }
    free(xdip);
    free(**itmp);
    free(*itmp);
    free(itmp);
    free(**xtmp);
    free(*xtmp);
    free(xtmp);
    free(**itmp2);
    free(*itmp2);
    free(itmp2);
//     for (i3=0; i3 < n3; i3++) {
// 	free(xtmp2[i3]);
//     }
//     free(xtmp2);
}

void pwsmooth3_lop(bool adj, bool add, 
		  int nin, int nout, float* trace, float *smooth)
/*< linear operator >*/
{
    int i1, i2, i3;

    sf_adjnull(adj,add,nin,nout,trace,smooth);

    if (adj) {
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    xtmp[i2][i3][i1] = smooth[i1+n1*(i2+n2*i3)];
		}
	    }
	}

	/* crossline */
	pwsmooth_init(ns2,n1,n3,order2,eps);
	for (i2=0; i2 < n2; i2++) {
	    pwsmooth_set(xdip[i2]);
	    pwsmooth_lop(true,false,n13,n13,itmp2[i2][0],xtmp[i2][0]);
	}
	pwsmooth_close();
	/* transpose */
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    itmp[i3][i2][i1] = itmp2[i2][i3][i1];
		}
	    }
	}
	/* inline */
	pwsmooth_init(ns1,n1,n2,order1,eps);
	for (i3=0; i3 < n3; i3++) {
	    pwsmooth_set(idip[i3]);
	    pwsmooth_lop(true,true,n12,n12,trace+i3*n12,itmp[i3][0]);
	}
	pwsmooth_close();
    } else {
	/* inline */
	pwsmooth_init(ns1,n1,n2,order1,eps);
	for (i3=0; i3 < n3; i3++) {
	    pwsmooth_set(idip[i3]);
	    pwsmooth_lop(false,false,n12,n12,trace+i3*n12,itmp[i3][0]);
	}
	pwsmooth_close();
	/* transpose */
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    itmp2[i2][i3][i1] = itmp[i3][i2][i1];
		}
	    }
	}
	/* crossline */
	pwsmooth_init(ns2,n1,n3,order2,eps);
	for (i2=0; i2 < n2; i2++) {
	    pwsmooth_set(xdip[i2]);
	    pwsmooth_lop(false,false,n13,n13,itmp2[i2][0],xtmp[i2][0]);
	}
	pwsmooth_close();
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    smooth[i1+n1*(i2+n2*i3)] += xtmp[i3][i2][i1];
		}
	    }
	}
    }
}
