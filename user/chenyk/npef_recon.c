/* 3-D missing data interpolation with nonstationary filters */
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

#include "npef_recon.h"

static float* filt;
static int nf1, nf2, nf3, nf4, nf5, n1, n2, n3, n4, n5;

void mmmult5_init (float* bb, 
		  int nf1t, int nf2t, int nf3t, int nf4t, int nf5t, 
		  int n1t,  int n2t,  int n3t,  int n4t,  int n5t) 
/*< initialize with a pointer to a matrix >*/
{
    filt = bb;
    nf1    = nf1t;
    nf2    = nf2t;
    nf3    = nf3t;
    nf4    = nf4t;
    nf5    = nf5t;
    n1     = n1t;
    n2     = n2t;
    n3     = n3t;
    n4     = n4t;
    n5     = n5t;
}

void mmmult5_lop (bool adj, bool add, 
		  int nx, int ny, float* mm, float* dd) 
/*< linear operator >*/
{
    int i1,i2,i3,i4,i5,j1,j2,j3,j4,j5;

    sf_adjnull (adj,add,nx,ny,mm,dd);

	    for (j5=0; j5 < n5; j5++) {
	    for (j4=0; j4 < n4; j4++) {  	
	    for (j3=0; j3 < n3; j3++) {
		for (j2=0; j2 < n2; j2++) {
		for (j1=0; j1 < n1; j1++) {
			    	for (i5=0; i5 < nf5; i5++) {
			    	for (i4=0; i4 < nf4; i4++) {			
			    	for (i3=0; i3 < nf3; i3++) {
					for (i2=0; i2 < nf2; i2++) {
				    for (i1=-nf1/2; i1 < (nf1+1)/2; i1++) {
			    	/* zero value boundary conditions */
				    if(j1+i1 <0 || j1+i1 >=n1 ||
				       j2+i2 <0 || j2+i2 >=n2 ||
				       j3+i3 <0 || j3+i3 >=n3 ||
				       j4+i4 <0 || j4+i4 >=n4 ||
				       j5+i5 <0 || j5+i5 >=n5)
				    	{continue;}
				    	if (adj)
				    	{
					  	mm[(j5+i5)*n4*n3*n2*n1+
					  	   (j4+i4)*n3*n2*n1+
					  	   (j3+i3)*n2*n1+
					  	   (j2+i2)*n1+
					  	   j1+i1]	
					  	+=filt[j5*n4*n3*n2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j4*n3*n2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j3*n2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j1*nf5*nf4*nf3*nf2*nf1+
					  	   i5*nf4*nf3*nf2*nf1+
					  	   i4*nf3*nf2*nf1+
					  	   i3*nf2*nf1+
					  	   i2*nf1+
					  	   i1+nf1/2]*
					  	dd[j5*n4*n3*n2*n1+
					  	   j4*n3*n2*n1+
					  	   j3*n2*n1+
					  	   j2*n1+
					  	   j1];
					  	}else{
					  	dd[j5*n4*n3*n2*n1+
					  	   j4*n3*n2*n1+
					  	   j3*n2*n1+
					  	   j2*n1+
					  	   j1]
					  	+=filt[j5*n4*n3*n2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j4*n3*n2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j3*n2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j2*n1*nf5*nf4*nf3*nf2*nf1+
					  	   j1*nf5*nf4*nf3*nf2*nf1+
					  	   i5*nf4*nf3*nf2*nf1+
					  	   i4*nf3*nf2*nf1+
					  	   i3*nf2*nf1+
					  	   i2*nf1+
					  	   i1+nf1/2]*
					  	mm[(j5+i5)*n4*n3*n2*n1+
					  	   (j4+i4)*n3*n2*n1+
					  	   (j3+i3)*n2*n1+
					  	   (j2+i2)*n1+
					  	   j1+i1];				  	
					  	
					  	
					  	}
				}
				}
			    }
				}
				}
			
		    }
		}
	    } 
	}
	}

}

void mmmult_close () 
/*< free filter memory >*/
{
    free (filt);
}


void nmis5(int niter         /* number of iterations */, 
	   int nf1, int nf2, int nf3, int nf4, int nf5,  /*NPEF filter size*/
	   int n1,  int n2,  int n3,  int n4, int n5,	 /*Data size*/
	   float *filt,
	   float *xx         /* model */, 
	   const bool *known /* mask for known data */,
	   float eps         /* regularization parameter */,
	   bool verb         /* verbosity flag */) 
/*< 3-D interpolation >*/
{
    int ix;
    float *dd;

    /*  regularized */
    dd = sf_floatalloc(n1*n2*n3*n4*n5);
    for (ix=0; ix < n1*n2*n3*n4*n5; ix++) {
	dd[ix]=0.;
    }
    
    mmmult5_init(filt, nf1, nf2, nf3, nf4, nf5, n1, n2, n3, n4, n5);
    sf_solver (mmmult5_lop, sf_cgstep, n1*n2*n3*n4*n5, n1*n2*n3*n4*n5, 
	       xx, dd, niter, "known", known, "x0", xx, "verb", verb, "end");

   // mmmult5_lop (0, 0, n1*n2*n3*n4*n5, n1*n2*n3*n4*n5, xx, dd);

 /*   for (ix=0; ix < n1*n2*n3*n4*n5; ix++) {
	xx[ix]=dd[ix];
    }*/

    free(dd);
    sf_cgstep_close();
    
}


