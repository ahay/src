/* 2-D Gaussian frequency-domain smoothing */
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

#include <math.h>

#include <rsf.h>

#include "flatpifilt.h"
#include "freqfilt4pi.h"
#include "Faddeeva.h"

//static float **shape;
static sf_complex **shape;

void flatpifilt_init(int n1, int n2     /* data size */, 
		 float d1, float d2 /* sampling */,
		 float passthr /* pass threshold */, 
		 float v_1 /* left cut velocity */,
		 float v_2 /* left full pass velocity */,
		 float v_3 /* right full pass velocity */,
		 float v_4 /* right cut velocity  */,
		 float eps)
/*< initialize (call freqfilt2 afterwards) >*/
{
    int ik, iw, nfft, nw;
    float dw, w, dk, k, k0, v_0, v_a, v_b, beta_1, beta_2, beta;
    double complex z, zl, zr, zc; // coefficients for erfi calculation
	
	//float beta=10.0; // bias coefficient
	//float v_a=0.0, v_b=2.0, v_0=1.0; // starting, ending, bias velocities
	//float eps = 0.001; // damp 

    /* determine frequency sampling (for real to complex FFT) */
    nfft = n1;
    if (n1%2) nfft++;
    nw = nfft/2+1;
    //dw = 2.*SF_PI/nfft;
    
    dw = 1./(n1*d1);
	//w0 = -0.5/d1;
    
    /* determine wavenumber sampling (for complex FFT) */
    //dk = 2.*SF_PI/n2;
    //k0 = -SF_PI;
    
    dk = 1./(n2*d2);
    k0 = -0.5/d2;

    // I propose to estimate beta 
    // based on four velocity points
    if (passthr == 999999999.999) {
	passthr = 1.0/0.001;
	}

    // calculating left slope
    beta_1 = (float)log(1.0/passthr);
    beta_1 /= (float)pow((v_2 - v_1),2);
    sf_warning("left beta=%f",beta_1);

    if(beta_1 > 30.0)
    {
	sf_warning("left beta > 30.0 - might crash (working on it)");
    }

    // calculating right slope
    beta_2 = (float)log(1.0/passthr);
    beta_2 /= (float)pow((v_3 - v_4),2);
    sf_warning("right beta=%f",beta_2);

    if(beta_2 > 30.0)
    {
	sf_warning("right beta > 30.0 - might crash (working on it)");
    }
    
    //shape = sf_floatalloc2(n2,nw);
	
	// need to make it double complex
	shape = sf_complexalloc2 (n2, nw);
 
    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	
	for (ik=0; ik < n2; ik++) {
	    k = k0+ik*dk;
	    
	    //Path-Integral Analytical Evaluation

	    // left slope
	    v_0 = v_2;
	    v_a = v_1;
	    v_b = v_2;
            beta = beta_1;
            zl = computepi(k,w,eps,v_0,v_a,v_b,beta);
	    
	    /*// check for NaN
            if ( creal(zl) != creal(zl)){
		if(!ch){
	    	sf_warning("Re(zl) = %f",creal(zl));
		ch = 1;
		}
	    }*/

            // right slope
            v_0 = v_3;
	    v_a = v_3;
	    v_b = v_4;
            beta = beta_2;
	    zr = computepi(k,w,eps,v_0,v_a,v_b,beta);
	    
            /*// check for NaN
            if ( creal(zr) != creal(zr)){
                if(!ch){
	    	sf_warning("Re(zr) = %f",creal(zr));
		ch = 1;
		}
	    }*/

            // center no weighting
            v_0 = 0.01;// any value - beta is zero
	    v_a = v_2;
	    v_b = v_3;
            beta = 0.0;
            zc = computepi(k,w,eps,v_0,v_a,v_b,beta);
	    
            /*// check for NaN
	    if ( creal(zc) != creal(zc)){
		if(!ch){
	    	sf_warning("Re(zc) = %f",creal(zc));
		ch = 1;
		}
	    }*/

            // sum
            z = zl + zc + zr;

	    //sf_warning("Re(z) = %f and Im(z) = %f",creal(z),cimag(z));
	    
	    // right now we assume that
	    // shape is double complex as well
	    // sf_cmplx ===== double complex
	    shape[iw][ik] = sf_cmplx(creal(z),cimag(z));
	    
	    //shape[iw][ik] = sf_cmplx(1.0,0.0);
	    
	}
    }

    freqfilt4pi_init(n1,n2,nw);
    freqfilt4pi_set(shape);
}

double complex computepi(float k, float w, float eps, float v_0, float v_a, float v_b, float beta)
/*< not for external use >*/
{ 

    double complex alpha, root, u_b, u_a, temp1, temp2, temp3, coeff, z; // coefficients for erfi calculation

    //computing coefficients for erfi
    alpha = (-1)*(k+eps)*(k+eps)*2.0*SF_PI/(16*(w+eps));
			
    root = csqrt(I*alpha - beta);

    //erfi arguments for v_a and v_b 
    u_b = v_b*root + v_0*beta/root;
    u_a = v_a*root + v_0*beta/root;
			
    //integral coefficient	
    coeff = cexp(-beta*v_0*v_0)*cexp(-beta*beta*v_0*v_0/(root*root))/root;
			
    temp1 = Faddeeva_Dawson(u_a,0);
    temp2 = Faddeeva_Dawson(u_b,0);
    
    temp3 = cexp(u_a*u_a);
			
    z = coeff*temp3*((cexp(u_b*u_b)*temp2/temp3) - temp1);
    z *= 2.0/(sqrt(SF_PI));

    return z;
}

void flatpifilt_close(void)
/*< free allocated storage >*/
{
    free(shape[0]);
    free(shape);
    freqfilt4pi_close();
}

/* 	$Id: gauss2.c 1131 2005-04-20 18:19:10Z fomels $	 */
