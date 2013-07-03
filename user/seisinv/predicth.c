/* Trace prediction in t-x plane with plane-wave destruction and average alonge offset axis*/
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

#include "predicth.h"
#include "pwd.h"

static int n1, n2, n3, nb, k2;
static sf_bands slv;
static float *diag, **offd, eps, eps2, **dip, *tt;
static pwd w1, w2;

static void stepper(bool adj /* adjoint flag */,
		    int i2   /* trace number */);

void predicth_init (int nx, int ny, int nh /* data size */, 
		   float e        /* regularization parameter */,
		   int nw         /* accuracy order */,
		   int k          /* radius */,
		   bool two       /* if two predictions */)
/*< initialize >*/
{
    n1 = nx;
    n2 = ny;
    n3 = nh;
    nb = 2*nw;

    eps = e;
    eps2 = e;

    slv = sf_banded_init (n1, nb);
    diag = sf_floatalloc (n1);
    offd = sf_floatalloc2 (n1,nb);

    w1 = pwd_init (n1, nw);
    w2 = two? pwd_init(n1,nw): NULL;

    tt = NULL;

    k2 = k;
    if (k2 > n2-1) sf_error("%s: k2=%d > n2-1=%d",__FILE__,k2,n2-1);
}

void predict_close (void)
/*< free allocated storage >*/
{
    sf_banded_close (slv);
    free (diag);
    free (*offd);
    free (offd);
    pwd_close (w1);
    if (NULL != w2) pwd_close (w2);
    if (NULL != tt) free(tt);
}

static void regularization(void)
/* fill diag and offd using regularization */ 
{
    int i1, ib;

    for (i1=0; i1 < n1; i1++) {
	diag[i1] = 6.*eps;
    	offd[0][i1] = -4.*eps;
    	offd[1][i1] = eps;
	for (ib=2; ib < nb; ib++) {
	    offd[ib][i1] = 0.0;
	}
    }

    diag[0] = diag[n1-1] = eps2+eps;
    diag[1] = diag[n1-2] = eps2+5.*eps;
    offd[0][0] = offd[0][n1-2] = -2.*eps;
}

void predict_step(bool adj            /* adjoint flag */,
		  bool forw           /* forward or backward */, 
		  float* trace        /* input/output trace */,
		  const float* pp    /* slope */)
/*< prediction step >*/
{
    float t0, t1, t2, t3;
    
    regularization();
    pwd_define (forw, w1, pp, diag, offd);
    sf_banded_define (slv, diag, offd);

    if (adj) sf_banded_solve (slv, trace);

    t0 = trace[0];
    t1 = trace[1];
    t2 = trace[n1-2];
    t3 = trace[n1-1];

    pwd_set (adj, w1, trace, trace, diag);

    trace[0] += eps2*t0;
    trace[1] += eps2*t1;
    trace[n1-2] += eps2*t2;
    trace[n1-1] += eps2*t3;

    if (!adj) sf_banded_solve (slv, trace);
}

void predict1_step(bool forw      /* forward or backward */, 
		   float* trace1  /* input trace */,
		   const float* pp /* slope */,
		   float* trace /* output trace */)
/*< prediction step from one trace >*/
{
    float t0, t1, t2, t3;
    
    regularization();

    pwd_define (forw, w1, pp, diag, offd);
    sf_banded_define (slv, diag, offd);

    t0 = trace1[0];
    t1 = trace1[1];
    t2 = trace1[n1-2];
    t3 = trace1[n1-1];

    pwd_set (false, w1, trace1, trace, diag);

    trace[0] += eps2*t0;
    trace[1] += eps2*t1;
    trace[n1-2] += eps2*t2;
    trace[n1-1] += eps2*t3;

    sf_banded_solve (slv, trace);
}

void predict2_step(bool forw1        /* forward or backward */, 
		   bool forw2,
		   float* trace1     /* input trace */,
		   float* trace2,
		   const float* pp1  /* slope */,
		   const float* pp2,
		   float *trace      /* output trace */)
/*< prediction step from two traces>*/
{
    int i1;
    float t0, t1, t2, t3;
    
    regularization();

    pwd_define (forw1, w1, pp1, diag, offd);
    pwd_define (forw2, w2, pp2, diag, offd);

    sf_banded_define (slv, diag, offd);

    t0 = 0.5*(trace1[0]+trace2[0]);
    t1 = 0.5*(trace1[1]+trace2[1]);
    t2 = 0.5*(trace1[n1-2]+trace2[n1-2]);
    t3 = 0.5*(trace1[n1-1]+trace2[n1-1]);

    pwd_set (false, w1, trace1, offd[0], trace);
    pwd_set (false, w2, trace2, offd[1], trace);

    for (i1=0; i1 < n1; i1++) {
	trace[i1] = offd[0][i1]+offd[1][i1];
    }

    trace[0] += eps2*t0;
    trace[1] += eps2*t1;
    trace[n1-2] += eps2*t2;
    trace[n1-1] += eps2*t3;

    sf_banded_solve (slv, trace);
}

void predict_set(float **dip1 /* dip field [n2][n1] */)
/*< set the local slopes for applying the linear operator >*/
{
    dip=dip1;
    if (NULL == tt) tt = sf_floatalloc(n1);
}

static void stepper(bool adj /* adjoint flag */,
		    int i2   /* trace number */)
{
    if (i2 < k2) {
	predict_step(adj,false,tt,dip[k2-1-i2]);
    } else if (i2 < n2+k2-1) {
	predict_step(adj,true,tt,dip[i2-k2]);
    } else {
	predict_step(adj,false,tt,dip[2*n2+k2-3-i2]);
    }
}

void predicth_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int i1, i2, i3, i;
    float t;

    if (nx != ny || nx != n1*n2*n3) sf_error("%s: Wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    sf_warning("n1=%d,n2=%d,n3=%d",n1,n2,n3);
    if (adj) {
        for (i3=0; i3 < n3; i3++) {
            for (i1=0; i1 < n1; i1++) {
	        tt[i1] = 0.;
            }

	    for (i2=n2-1; i2 >= 0; i2--) {
	        predict_step(true,true,tt,dip[i2]);
	        for (i1=0; i1 < n1; i1++) {
		    tt[i1] += yy[i1+i2*n1+i3*n1*n2];
	        }
	        for (i1=0; i1 < n1; i1++) {
		    xx[i1+i2*n1+i3*n1*n2] += tt[i1];
	        }
             }
	}
        for (i2=0; i2 < n2; i2++) {
            for (i1=0; i1 < n1; i1++) {
                t=0.;
                for (i3=n3-1,i=1; i3 >= 0; i3--,i++) {
                    t += xx[i3*n2*n1+i2*n1+i1];
                    xx[i3*n2*n1+i2*n1+i1] = t/i;
                }
            }
         }
    } else {

       for (i2=0; i2 < n2; i2++) {
           for (i1=0; i1 < n1; i1++) {
               t=0.;
               for (i3=0,i=1; i3 <= n3-1; i3++,i++) {
                   t += xx[i3*n2*n1+i2*n1+i1];
                   xx[i3*n2*n1+i2*n1+i1] = t/i;
               }
           }
        }

        for (i3=0; i3 < n3; i3++) {
            for (i1=0; i1 < n1; i1++) {
	        tt[i1] = 0.;
            }

	    for (i2=0; i2 < n2; i2++) {
	        for (i1=0; i1 < n1; i1++) {
		    tt[i1] += xx[i1+i2*n1+i3*n1*n2];
	        }
	        for (i1=0; i1 < n1; i1++) {
		    yy[i1+i2*n1+i3*n1*n2] += tt[i1];
	        }
	        predict_step(false,true,tt,dip[i2]);
	    }
         }
    }
}

void predicter_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int i1, i2;

    if (nx != ny || nx != n1*(n2+2*k2)) 
	sf_error("%s: Wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (i1=0; i1 < n1; i1++) {
	tt[i1] = 0.;
    }

    if (adj) {
	for (i2=n2+2*k2-1; i2 >= 0; i2--) {
	    stepper(true,i2);
	    for (i1=0; i1 < n1; i1++) {
		tt[i1] += yy[i1+i2*n1];
	    }
	    for (i1=0; i1 < n1; i1++) {
		xx[i1+i2*n1] += tt[i1];
	    }
	}
    } else {
	for (i2=0; i2 < n2+2*k2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		tt[i1] += xx[i1+i2*n1];
	    }
	    for (i1=0; i1 < n1; i1++) {
		yy[i1+i2*n1] += tt[i1];
	    }
	    stepper(false,i2);
	}
    }
}

void subtracter_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int i1, i2, j2, m2;

    if (nx != ny || nx != n1*(n2+2*k2)) 
	sf_error("%s: Wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    if (adj) {
	for (j2=0; j2 < n2+2*k2; j2++) {
	    i2=j2+k2;
	    if (i2 < n2+2*k2) {
		for (i1=0; i1 < n1; i1++) {
		    tt[i1] = yy[i1+i2*n1];
		}
		for (m2=i2-1; m2 >= j2; m2--) {
		    stepper(true,m2);
		}
		for (i1=0; i1 < n1; i1++) {
		    xx[i1+j2*n1] += yy[i1+j2*n1]-tt[i1];
		}
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    xx[i1+j2*n1] += yy[i1+j2*n1];
		}
	    }
	}
    } else {
	for (i2=0; i2 < n2+2*k2; i2++) { 
	    j2=i2-k2;
	    if (j2 >=0) {
		for (i1=0; i1 < n1; i1++) {
		    tt[i1] = xx[i1+j2*n1];
		}
		for (m2=j2; m2 < i2; m2++) {
		    stepper(false,m2);
		}
		for (i1=0; i1 < n1; i1++) {
		    yy[i1+i2*n1] += xx[i1+i2*n1]-tt[i1];
		}
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    yy[i1+i2*n1] += xx[i1+i2*n1];
		}
	    }
	}
    }
}

void subtract_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int i1, i2, j2, m2;

    if (nx != ny || nx != n1*n2) sf_error("%s: Wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    if (adj) {
	for (j2=0; j2 < n2; j2++) {
	    i2=j2+k2;
	    if (i2 < n2) {
		for (i1=0; i1 < n1; i1++) {
		    tt[i1] = yy[i1+i2*n1];
		}
		for (m2=i2-1; m2 >= j2; m2--) {
		    predict_step(true,true,tt,dip[m2]);
		}
		for (i1=0; i1 < n1; i1++) {
		    xx[i1+j2*n1] += yy[i1+j2*n1]-tt[i1];
		}
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    xx[i1+j2*n1] += yy[i1+j2*n1];
		}
	    }
	}
    } else {
	for (i2=0; i2 < n2; i2++) { 
	    j2=i2-k2;
	    if (j2 >=0) {
		for (i1=0; i1 < n1; i1++) {
		    tt[i1] = xx[i1+j2*n1];
		}
		for (m2=j2; m2 < i2; m2++) {
		    predict_step(false,true,tt,dip[m2]);
		}
		for (i1=0; i1 < n1; i1++) {
		    yy[i1+i2*n1] += xx[i1+i2*n1]-tt[i1];
		}
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    yy[i1+i2*n1] += xx[i1+i2*n1];
		}
	    }
	}
    }
}



