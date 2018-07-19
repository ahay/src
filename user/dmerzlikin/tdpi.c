/* 3-D Analytical Path-Summation Integral Evaluation */
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
/*^*/

#include "Faddeeva.h"
#include "fft3.h"

static sf_complex *c;
static float *f, beta_1, beta_2, v_1, v_2, v_3, v_4, eps1;
static int on1, on2, on3, nn1, nn2, nn3, nk;
static float dw, dkx,dky, k0x,k0y;

double complex computepi(float kx, float ky, float w, float eps, float v_0, float v_a, float v_b, float beta)
/*< not for external use >*/
{ 

    double complex alpha, root, u_b, u_a, temp1, temp2, temp3, coeff, z; // coefficients for erfi calculation
    int ch=0;

    //Path-Integral Analytical Evaluation Method 1

    /*//computing coefficients for erfi
    alpha = (-1)*((kx*kx+eps) + (ky*ky+eps))*2.0*SF_PI/(16*(w+eps));
			
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
    z *= 2.0/(sqrt(SF_PI));*/

    //Path-Integral Analytical Evaluation Method 2
	    
    //computing coefficients for erfi
    kx *= 2.*SF_PI; kx *= kx;
    ky *= 2.*SF_PI; ky *= ky;
    w *= 2.*SF_PI;

    alpha = (-1)*((kx+eps) + (ky+eps))/(16*(w+eps));
			
    root = csqrt(I*alpha - beta);

    //erfi arguments for v_a and v_b 
    u_b = v_b*root + v_0*beta/root;
    u_a = v_a*root + v_0*beta/root;
			
    //integral coefficient	
    coeff = cexp(-beta*v_0*v_0)*cexp(-beta*beta*v_0*v_0/(root*root))/root;

    if(creal(coeff) != creal(coeff)){
	if(!ch){
		sf_warning("tdpi.c: unstable coefficient");
		ch = 1;
	}		
    }

    temp1 = coeff*Faddeeva_erfi(u_a,0);
    temp2 = coeff*Faddeeva_erfi(u_b,0);

    if(creal(temp1) != creal(temp1)){
	if(!ch){
		sf_warning("tdpi.c: unstable erfi");
		ch = 1;
	}		
    }

    if(creal(temp2) != creal(temp2)){
	if(!ch){
		sf_warning("tdpi.c: unstable erfi");
		ch = 1;
	}		
    }


    z = temp2 - temp1;

    return z;

}

void tdpi_init(int n1, int n2, int n3  /* data size */, 
		 float d1, float d2, float d3 /* sampling */,
		 float passthr /* pass threshold */, 
		 float v1 /* left cut velocity */,
		 float v2 /* left full pass velocity */,
		 float v3 /* right full pass velocity */,
		 float v4 /* right cut velocity  */,
		 float eps)
/*< initialize path-summation integral and ffts >*/
{

    int nw, ch=0;

    on1 = n1;
    on2 = n2;
    on3 = n3;

    //if ( (n1 != nn1) || (n2 != nn2) || (n3 != nn3)) sf_warning("we will need to account for padding.");

    nk = fft3_init(false,1,on1,on2,on3,&nn1,&nn2,&nn3);
    
    c = sf_complexalloc(nk*nn2*nn3);
    f = sf_floatalloc(nn1*nn2*nn3);

    /* size should have been computed in fft3 init */
    nw = nk;
    dw = 1./(nn1*d1);
    dkx = 1./(nn2*d2);
    k0x = -0.5/d2;
    dky = 1./(nn3*d3);
    k0y = -0.5/d3;

    v_1 = v1;
    v_2 = v2;
    v_3 = v3;
    v_4 = v4;

    eps1 = eps;

    // I propose to estimate beta 
    // based on four velocity points
    if (passthr == 999999999.999) {
	passthr = 1.0/0.001;
	}

    // calculating left slope
    beta_1 = (float)log(1.0/passthr);
    beta_1 /= (float)pow((v_2 - v_1),2);
    //sf_warning("left beta=%f",beta_1);

    //if(beta_1 > 30.0)
    //{
	//sf_warning("left beta > 30.0 - might crash (working on it)");
    //}

    // calculating right slope
    beta_2 = (float)log(1.0/passthr);
    beta_2 /= (float)pow((v_3 - v_4),2);
    //sf_warning("right beta=%f",beta_2);

    //if(beta_2 > 30.0)
    //{
	//sf_warning("right beta > 30.0 - might crash (working on it)");
    //}

}

void tdpi_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< apply path-summation integral and ffts >*/
{

    int ix, iy, iz, index, ch=0;
    sf_file fftFile;
    bool fftDo;
    float v_0, v_a, v_b, beta, kx, ky, w;
    sf_complex temp;
    double complex z, zl, zr, zc; // coefficients for erfi calculation
    float complex zf;

    sf_adjnull(adj,add,nx,ny,x,y);

    /* zeroing intermediate allocations */
    for(ix=0; ix<nk*nn2*nn3; ix++){

	c[ix]=0.0;

    }

    for(ix=0; ix<nn1*nn2*nn3; ix++){

	f[ix]=0.0;

    }
    
    if (!sf_getbool("fftDo",&fftDo)) fftDo=false;

    if (fftDo) {

    	fftFile = sf_output("fft");
        sf_putint(fftFile,"n1",nk);
	sf_putint(fftFile,"n2",nn2);
	sf_putint(fftFile,"n3",nn3);
        sf_putfloat(fftFile,"o1",0.0);
	sf_putfloat(fftFile,"o2",k0x);
        sf_putfloat(fftFile,"o3",k0y);
        sf_putfloat(fftFile,"d1",dw);
	sf_putfloat(fftFile,"d2",dkx);
        sf_putfloat(fftFile,"d3",dky);
        sf_settype (fftFile,SF_COMPLEX);

    } else fftFile=NULL;

    for (ix=0; ix < nn1*nn2*nn3; ix++) {
	    f[ix]=0.;
    }

    for (iy=0; iy < on3; iy++) {
	for (ix=0; ix < on2; ix++) {
		for (iz=0; iz < on1; iz++) {
			f[(iy*nn2+ix)*nn1 + iz] = adj? y[(iy*on2+ix)*on1 + iz]: x[(iy*on2+ix)*on1 + iz];
		}
    	}

    }

    fft3(f,c);

    for (iy=0; iy < nn3; iy++) {

	ky = k0y+iy*dky;

	for (ix=0; ix < nn2; ix++) {

		kx = k0x+ix*dkx;

		for (iz=0; iz < nk; iz++) {

			w = iz*dw;

	                index = (iy*nn2+ix)*nk + iz;

			/* Path-Integral Analytical Evaluation */

		    	/* left slope */
		    	v_0 = v_2;
		    	v_a = v_1;
		    	v_b = v_2;
	            	beta = beta_1;
	            	zl = computepi(kx,ky,w,eps1,v_0,v_a,v_b,beta);
	    
		    	/* check for NaN in left path-summation integral */
	            	if ( creal(zl) != creal(zl)){
			
				if(!ch){
		    			sf_warning("tdpi.c: left p-si is unstable: try different range");
					ch = 1;
				}

		    	}

	            	/* right slope */
	            	v_0 = v_3;
		    	v_a = v_3;
		    	v_b = v_4;
	            	beta = beta_2;
		    	zr = computepi(kx,ky,w,eps1,v_0,v_a,v_b,beta);
	    
	            	/* check for NaN in right path-summation integral */
	            	if ( creal(zr) != creal(zr)){
			
				if(!ch){
		    			sf_warning("tdpi.c: right p-si is unstable: try different range");
					ch = 1;
				}

		    	}

	            	// center no weighting
	            	v_0 = 0.01;// any value - beta is zero
		    	v_a = v_2;
		    	v_b = v_3;
	            	beta = 0.0;
	            	zc = computepi(kx,ky,w,eps1,v_0,v_a,v_b,beta);
	    
	            	/* check for NaN in middle path-summation integral */
	            	if ( creal(zc) != creal(zc)){
			
				if(!ch){
		    			sf_warning("tdpi.c: middle p-si is unstable: try different range");
					ch = 1;
				}

		    	}
	
	            	/* combine all the parts */
	            	z = zl + zc + zr;

			zf = (float complex) z;

			/* check for NaN in middle path-summation integral */
	            	if ( crealf(zf) != crealf(zf)){
			
				if(!ch){
		    			sf_warning("tdpi.c: combination of p-si is unstable: try different range");
					ch = 1;
				}

		    	}

		if (adj){
		    
                    temp = sf_cmplx(crealf(z),(-1.0)*cimagf(z));

		} else {
		    
                    temp = sf_cmplx(crealf(z),cimagf(z));

		}

                c[index] = c[index]*temp;

		}
    	}
    }

    if(fftDo)
	{
	for (iy=0; iy < nn3; iy++){
		for (ix=0; ix < nn2; ix++){
			sf_complexwrite(c+(iy*nn2+ix)*nk,nk,fftFile);
		}
	}
    }

    ifft3_allocate(c);

    ifft3(f,c);

    for (iy=0; iy < on3; iy++) {
	for (ix=0; ix < on2; ix++) {
		for (iz=0; iz < on1; iz++) {
			if (adj == false){
			
				y[(iy*on2+ix)*on1 + iz] += f[(iy*nn2+ix)*nn1 + iz];

			} else {

				x[(iy*on2+ix)*on1 + iz] += f[(iy*nn2+ix)*nn1 + iz];

			}
		}
    	}

    }

}

void tdpi_close(void)
/*< free allocated storage >*/
{
    free(f);
    free(c);
    sf_warning("tdpi: memory free");
}
