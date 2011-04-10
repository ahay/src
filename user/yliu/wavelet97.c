/* Digital 9/7 biorthogonal wavelet transform */
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

#include "wavelet97.h"

static int nt;
static float *t;
static bool inv;

static void biorthogonal(bool adj) 
/* Lifting CDF 9/7 biorthogonal wavelet transform in place */
{
    int i,j;
    float a;

    if (adj) {
	for (j=nt/2; j >= 1; j /= 2) {
	    if (inv) {                            /*reverse dwt9/7 transform*/
                a= 1.230174105f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]  /= a;
	        }
	        t[0] /= a;                         /*left boundary*/
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] *= a;
		    /* Undo Scale */
	        }
	        if (i+j < nt) t[i+j] *= a;         /*right boundary*/  

                a= -0.4435068522f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]   += (t[i+j]+t[i-j])*a;
		    /* Undo Update 2 */
	        }
	        t[0] += 2*a*t[j];                  /*left boundary*/

	        a = -0.8829110762f;
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] += (t[i]+t[i+2*j])*a;
		    /* Undo Predict 2 */
	        }	 
	        if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
                    /* Undo Step 2 */

                a= 0.05298011854f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]   += (t[i+j]+t[i-j])*a;
		    /* Undo Update 1 */
	        }
	        t[0] += 2*a*t[j];                  /*left boundary*/ 

	        a = 1.586134342f;
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] += (t[i]+t[i+2*j])*a;
		    /* Undo Predict 1 */
	        }	 
	        if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
                    /* Undo Step 1 */


	    } else {                               /*adjoint transform*/
                a= 1.230174105f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i]  *= a;
	        }
	        t[0] *= a;                         /*left boundary*/
	        for (i=0; i < nt-2*j; i += 2*j) {
		    t[i+j] /= a;
		    /* Undo Scale */
	        }
	        if (i+j < nt) t[i+j] /= a;         /*right boundary*/  

                a= -0.4435068522f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i+j]   -= t[i]*a;
                    t[i-j]   -= t[i]*a;
		    /* Undo Update 2 */
	        }
	        t[j] -= 2*a*t[0];                  /*left boundary*/

	        a = -0.8829110762f;
	        for (i=0; i < nt-2*j; i += 2*j) {
                    t[i]     -= t[i+j]*a;
                    t[i+2*j] -= t[i+j]*a;
		    /* Undo Predict 2 */
	        }	 
	        if (i+j < nt) t[i] -= 2*a*t[i+j];  /*right boundary*/  
                    /* Undo Step 2 */

                a= 0.05298011854f;
	        for (i=2*j; i < nt-j; i += 2*j) {
		    t[i+j]   -= t[i]*a;
                    t[i-j]   -= t[i]*a;
		    /* Undo Update 1 */
	        }
	        t[j] -= 2*a*t[0];                  /*left boundary*/ 

	        a = 1.586134342f;
	        for (i=0; i < nt-2*j; i += 2*j) {
                    t[i]     -= t[i+j]*a;
                    t[i+2*j] -= t[i+j]*a;
		    /* Undo Predict 1 */
	        }	 
	        if (i+j < nt) t[i] -= 2*a*t[i+j];  /*right boundary*/  
                    /* Undo Step 1 */

	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {        /*different scale*/
	    a = -1.586134342f;
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] += (t[i]+t[i+2*j])*a;
		/* Predict 1 */
	    }	 
	    if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
 
            a= -0.05298011854f;
	    t[0] += 2*a*t[j];                  /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]   += (t[i+j]+t[i-j])*a;
		/* Update 1 */
	    }
                /* Step 1 */

	    a = 0.8829110762f;
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] += (t[i]+t[i+2*j])*a;
		/* Predict 2 */
	    }	 
	    if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
 
            a= 0.4435068522f;
	    t[0] += 2*a*t[j];                  /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]   += (t[i+j]+t[i-j])*a;
		/* Update 2 */
	    }
                /* Step 2 */

            a= 1/(1.230174105f);
	    for (i=0; i < nt-2*j; i += 2*j) {
		t[i+j] *= a;
	    }	 
	    if (i+j < nt) t[i+j] *= a;         /*right boundary*/  
	    t[0] /= a;                         /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
		t[i]  /= a;
		/* Scale */
	    }
	}
    }

}


void wavelet97_init(int n /* data size */, bool inv1) 
/*< allocate space >*/
{
    inv = inv1;
    for (nt=1; nt < n; nt *= 2) ;  
    /*get nt, nt is 2^n and less than n */
    t = sf_floatalloc(nt);

}

void wavelet97_close(void) 
/*< deallocate space >*/
{
    free (t);
}

void wavelet97_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int it, i, j;
    sf_adjnull (adj,add,nx,ny,x,y);

    if (adj) {
        t[0] = y[0];
        it = 1;
        for (j=nt/2;j>=1;j/=2) {
            for (i=0;i<nt-j;i+=2*j) {
                if (it < ny) {
                   t[i+j]=y[it];
                   it++;
                } else {
                    t[i+j]=0.;
                }
            }
        }
    } else {
        for (it=0; it < nx; it++) {
            t[it]=x[it];
        }
        for (it=nx; it < nt; it++) {
            t[it] = 0.;
        }
    }
    biorthogonal(adj);

    if (adj) {
        for (it=0;it <nx;it++) {
            x[it] +=t[it];
        }
    } else {
        y[0] += t[0];
        it = 1;
        for (j=nt/2; j >= 1; j /=2) {
            for (i=0; i < nt-j; i +=2*j) {
               y[it] +=t[i+j];
               it++;
               if( it>= ny) return;
            }
        }
    }

}
