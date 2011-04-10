/* Digital freqlet transform */
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

#include "freqlet97.h"

static int nt;
static bool inv;
static sf_complex *t, *z;

static void biorthogonal(bool adj) 
/* Lifting CDF 9/7 biorthogonal wavelet transform in place */
{
    int i, j;
    sf_complex z0;
    float a;

    if (adj) {

	for (j=nt/2; j >= 1; j /= 2) {
            z0 = z[j];
	    if (inv) {                            /*reverse dwt9/7 transform*/
                a= 1.230174105f;
	        for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i]  /= a;
#else
	            t[i] = sf_crmul(t[i],(1/a));
#endif
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[0] /= a;                         /*left boundary*/
#else
		t[0] = sf_crmul(t[0],(1/a));
#endif
	        for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] *= a;
#else
	            t[i+j] = sf_crmul(t[i+j],a);
#endif
		    /* Undo Scale */
	        }
#ifdef SF_HAS_COMPLEX_H
	        if (i+j < nt) t[i+j] *= a;         /*right boundary*/  
#else
	        if (i+j < nt) t[i+j] = sf_crmul(t[i+j],a);
#endif
                a= -0.4435068522f;
	        for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i]   += (t[i+j]/z0+t[i-j]*z0)*a;
#else
		    t[i] = sf_cadd(t[i],
				   sf_crmul(
				       sf_cadd(
					   sf_cmul(t[i+j],sf_conjf(z0)),
					   sf_cmul(t[i-j],z0)),
				       a));
#endif
		    /* Undo Update 2 */
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[0] += 2*a*t[j]/z0;                  /*left boundary*/
#else
		t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(sf_cmul(t[j],sf_conjf(z0)),a),2));
#endif
	        a = -0.8829110762f;
	        for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += (t[i]*z0+t[i+2*j]/z0)*a;
#else
		    t[i+j] = sf_cadd(t[i+j],
				   sf_crmul(
				       sf_cadd(
					   sf_cmul(t[i],z0),
					   sf_cmul(t[i+2*j],sf_conjf(z0))),
				       a));
#endif
		    /* Undo Predict 2 */
	        }
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < nt) t[i+j] += 2*a*t[i]*z0;  /*right boundary*/  
#else
		if (i+j < nt) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(sf_cmul(t[i],z0),a),2));
#endif
                    /* Undo Step 2 */

                a= 0.05298011854f;

	        for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i]   += (t[i+j]/z0+t[i-j]*z0)*a;
#else
		    t[i] = sf_cadd(t[i],
				   sf_crmul(
				       sf_cadd(
					   sf_cmul(t[i+j],sf_conjf(z0)),
					   sf_cmul(t[i-j],z0)),
				       a));
#endif
		    /* Undo Update 1 */
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[0] += 2*a*t[j]/z0;                  /*left boundary*/
#else
		t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(sf_cmul(t[j],sf_conjf(z0)),a),2));
#endif

	        a = 1.586134342f;
	        for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] += (t[i]*z0+t[i+2*j]/z0)*a;
#else
		    t[i+j] = sf_cadd(t[i+j],
				   sf_crmul(
				       sf_cadd(
					   sf_cmul(t[i],z0),
					   sf_cmul(t[i+2*j],sf_conjf(z0))),
				       a));
#endif
		    /* Undo Predict 1 */
	        }	 
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < nt) t[i+j] += 2*a*t[i]*z0;  /*right boundary*/  
#else
		if (i+j < nt) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(sf_cmul(t[i],z0),a),2));
#endif
                    /* Undo Step 1 */


	    } else {                               /*adjoint transform*/
                a= 1.230174105f;
	        for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i]  *= a;
#else
	            t[i] = sf_crmul(t[i],a);
#endif
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[0] *= a;                         /*left boundary*/
#else
		t[0] = sf_crmul(t[0],a);
#endif
	        for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j] /= a;
#else
	            t[i+j] = sf_crmul(t[i+j],(1/a));
#endif
		    /* Undo Scale */
	        }
#ifdef SF_HAS_COMPLEX_H
	        if (i+j < nt) t[i+j] /= a;         /*right boundary*/  
#else
	        if (i+j < nt) t[i+j] = sf_crmul(t[i+j],(1/a));
#endif
                a= -0.4435068522f;
	        for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j]   -= t[i]*z0*a;
                    t[i-j]   -= t[i]/z0*a;
#else
		    t[i+j] = sf_cadd(t[i+j],
				   sf_crmul(sf_cmul(t[i],z0),(-a)));
		    t[i-j] = sf_cadd(t[i-j],
				   sf_crmul(sf_cmul(t[i],sf_conjf(z0)),(-a)));
#endif
		    /* Undo Update 2 */
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[j] -= 2*a*t[0]*z0;                  /*left boundary*/
#else
		t[j] = sf_cadd(t[j],sf_crmul(sf_crmul(sf_cmul(t[0],z0),a),-2));
#endif
	        a = -0.8829110762f;
	        for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
                    t[i]     -= (t[i+j]/z0)*a;
                    t[i+2*j] -= (t[i+j]*z0)*a;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(sf_cmul(t[i+j],sf_conjf(z0)),(-a)));
		    t[i+2*j] = sf_cadd(t[i+2*j],sf_crmul(sf_cmul(t[i+j],z0),(-a)));
#endif
		    /* Undo Predict 2 */
	        }
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < nt) t[i] -= 2*a*t[i+j]/z0;  /*right boundary*/  
#else
		if (i+j < nt) t[i] = sf_cadd(t[i],sf_crmul(sf_crmul(sf_cmul(t[i],sf_conjf(z0)),a),-2));
#endif
                    /* Undo Step 2 */
 
                a= 0.05298011854f;
	        for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		    t[i+j]   -= t[i]*z0*a;
                    t[i-j]   -= t[i]/z0*a;
#else
		    t[i+j] = sf_cadd(t[i+j],
				   sf_crmul(sf_cmul(t[i],z0),(-a)));
		    t[i-j] = sf_cadd(t[i-j],
				   sf_crmul(sf_cmul(t[i],sf_conjf(z0)),(-a)));
#endif		    /* Undo Update 1 */
	        }
#ifdef SF_HAS_COMPLEX_H
	        t[j] -= 2*a*t[0]*z0;                  /*left boundary*/
#else
		t[j] = sf_cadd(t[j],sf_crmul(sf_crmul(sf_cmul(t[0],z0),a),-2));
#endif
	        a = 1.586134342f;
	        for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
                    t[i]     -= (t[i+j]/z0)*a;
                    t[i+2*j] -= (t[i+j]*z0)*a;
#else
		    t[i] = sf_cadd(t[i],sf_crmul(sf_cmul(t[i+j],sf_conjf(z0)),(-a)));
		    t[i+2*j] = sf_cadd(t[i+2*j],sf_crmul(sf_cmul(t[i+j],z0),(-a)));
#endif
		    /* Undo Predict 1 */
	        }	 
#ifdef SF_HAS_COMPLEX_H	 
	        if (i+j < nt) t[i] -= 2*a*t[i+j]/z0;  /*right boundary*/  
#else
		if (i+j < nt) t[i] = sf_cadd(t[i],sf_crmul(sf_crmul(sf_cmul(t[i],sf_conjf(z0)),a),-2));
#endif                    /* Undo Step 1 */
	    }
	}
    } else {
	for (j=1; j <= nt/2; j *= 2) {        /*different scale*/
	    z0=z[j];
	    a = -1.586134342f;
	    for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
	  	 t[i+j] += (t[i]*z0+t[i+2*j]/z0)*a;
#else
	   	 t[i+j] = sf_cadd(t[i+j],
		 	      sf_crmul(
		 	          sf_cadd(
		   	              sf_cmul(t[i],z0),
				      sf_cmul(t[i+2*j],sf_conjf(z0))),
				  a));
#endif
//	    a = -1.586134342f;
//	    for (i=0; i < nt-2*j; i += 2*j) {
//		t[i+j] += (t[i]+t[i+2*j])*a;
		/* Predict 1 */
	    }	 

#ifdef SF_HAS_COMPLEX_H	 
            if (i+j < nt) t[i+j] += 2*a*t[i]*z0;  /*right boundary*/  
#else
	    if (i+j < nt) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(sf_cmul(t[i],z0),a),2));
#endif
//	    if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
 
            a= -0.05298011854f;
#ifdef SF_HAS_COMPLEX_H
	    t[0] += 2*a*t[j]/z0;                  /*left boundary*/
#else
	    t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(sf_cmul(t[j],sf_conjf(z0)),a),2));
#endif
//	    t[0] += 2*a*t[j];                  /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i]   += (t[i+j]/z0+t[i-j]*z0)*a;
#else
		t[i] = sf_cadd(t[i],
			   sf_crmul(
			       sf_cadd(
				   sf_cmul(t[i+j],sf_conjf(z0)),
				   sf_cmul(t[i-j],z0)),
			       a));
#endif
//		t[i]   += (t[i+j]+t[i-j])*a;
		/* Update 1 */
	    }
                /* Step 1 */

	    a = 0.8829110762f;
	    for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
	  	 t[i+j] += (t[i]*z0+t[i+2*j]/z0)*a;
#else
	   	 t[i+j] = sf_cadd(t[i+j],
		 	      sf_crmul(
		 	          sf_cadd(
		   	              sf_cmul(t[i],z0),
				      sf_cmul(t[i+2*j],sf_conjf(z0))),
				  a));
#endif
//		t[i+j] += (t[i]+t[i+2*j])*a;
		/* Predict 2 */
	    }	 
#ifdef SF_HAS_COMPLEX_H	 
            if (i+j < nt) t[i+j] += 2*a*t[i]*z0;  /*right boundary*/  
#else
	    if (i+j < nt) t[i+j] = sf_cadd(t[i+j],sf_crmul(sf_crmul(sf_cmul(t[i],z0),a),2));
#endif
//	    if (i+j < nt) t[i+j] += 2*a*t[i];  /*right boundary*/  
 
            a= 0.4435068522f;
#ifdef SF_HAS_COMPLEX_H
	    t[0] += 2*a*t[j]/z0;                  /*left boundary*/
#else
	    t[0] = sf_cadd(t[0],sf_crmul(sf_crmul(sf_cmul(t[j],sf_conjf(z0)),a),2));
#endif
//	    t[0] += 2*a*t[j];                  /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i]   += (t[i+j]/z0+t[i-j]*z0)*a;
#else
		t[i] = sf_cadd(t[i],
			   sf_crmul(
			       sf_cadd(
				   sf_cmul(t[i+j],sf_conjf(z0)),
				   sf_cmul(t[i-j],z0)),
			       a));
#endif
//		t[i]   += (t[i+j]+t[i-j])*a;
		/* Update 2 */
	    }
                /* Step 2 */

            a= 1/(1.230174105f);
	    for (i=0; i < nt-2*j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i+j] *= a;
#else
	        t[i+j] = sf_crmul(t[i+j],a);
#endif
//		t[i+j] *= a;
	    }	 
#ifdef SF_HAS_COMPLEX_H
	    if (i+j < nt) t[i+j] *= a;         /*right boundary*/  
#else
	    if (i+j < nt) t[i+j] = sf_crmul(t[i+j],a);
#endif
//	    if (i+j < nt) t[i+j] *= a;         /*right boundary*/  
#ifdef SF_HAS_COMPLEX_H
	    t[0] /= a;                         /*left boundary*/
#else
	    t[0] = sf_crmul(t[0],(1/a));
#endif
//	    t[0] /= a;                         /*left boundary*/
	    for (i=2*j; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		t[i]  /= a;
#else
	        t[i] = sf_crmul(t[i],(1/a));
#endif
//		t[i]  /= a;
		/* Scale */
	    }
	}
    }
}

void freqlet97_init(int n /* data size */, bool inv1) 
/*< allocate space >*/
{
    int j;
    inv = inv1;

    for (nt=1; nt < n; nt *= 2) ;
    t = sf_complexalloc(nt);
    z = sf_complexalloc(nt);
    
    for (j=1; j <= nt/2; j *= 2) {
	z[j] = sf_cmplx(1.,0.); 
    }
}

void freqlet97_set(float w0)
/*< set frequency >*/
{
    int j;

    for (j=1; j <= nt/2; j *= 2) {
	z[j] = sf_cmplx(cosf(w0*j),sinf(w0*j));
    }
}

void freqlet97_close(void) 
/*< deallocate space >*/
{
    free (t);
    free (z);
}

void freqlet97_lop(bool adj, bool add, int nx, int ny, 
		 sf_complex *x, sf_complex *y)
/*< linear operator >*/
{
    int it, i, j;

    sf_cadjnull (adj,add,nx,ny,x,y);

    if (adj) {
	t[0] = y[0];
	it = 1;
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		if (it < ny) {
		    t[i+j]=y[it];
		    it++;
		} else {
		    t[i+j]=sf_cmplx(0.,0.);
		}
	    }	    	    
	}

    } else {
	for (it=0; it < nx; it++) {
	    t[it]=x[it];
	}
	for (it=nx; it < nt; it++) {
	    t[it] = sf_cmplx(0.,0.);
	}
    }

    biorthogonal(adj);    

    if (adj) {
	for (it=0; it < nx; it++) {
#ifdef SF_HAS_COMPLEX_H
	    x[it] += t[it];
#else
	    x[it] = sf_cadd(x[it],t[it]);
#endif
	}
    } else {

#ifdef SF_HAS_COMPLEX_H
	y[0] += t[0];
#else
	y[0] = sf_cadd(y[0],t[0]);
#endif
	it = 1;
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
#ifdef SF_HAS_COMPLEX_H
		y[it] += t[i+j];
#else
		y[it] = sf_cadd(y[it],t[i+j]);
#endif
		it++;
		if (it >= ny) return;
	    }	    	    
	}
    }
}



