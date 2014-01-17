/* 1D signal analysis by using Duffing differential equation solved by 4th order Runge-Kutta method. 
Duffing equation: x''/(omega^2)+0.5 x'/omega-x+x^3=gamma cos(omega t+phi)+kxi R(t)
*/
/*
  Copyright (C) 2014 Jilin University
  
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
#include <stdio.h>
#include <math.h>
float function1(float t,float x,float y,float omega)
{
    return omega*y;
}

float function2(float t,float x,float y,
		float Rd,float kxi,float gamma,
		float omega,int pow1,int pow2,float phi)
{
    float y1;
    y1=omega*(powf(x,pow1)-powf(x,pow2)-
	      0.5*y+kxi*Rd+gamma*cos(omega*t+phi*3.1415926));
    return y1;
}

float function3(float t,float x,float y,
		float Rd,float kxi,float gamma,
		float omega,int pow1,int pow2,float restore)
{
    float y1;
    y1=omega*(powf(x,pow1)-powf(x,pow2)-0.5*y+kxi*Rd+gamma*restore);
    return y1;
}

int main(int argc,char *argv[])
{
    float x0,y0,x00,y00,t0,d1,d2,d3,o1,o2,o3,k11,k21,k31,k41,k12,k22,k32,k42;
    float kxi,gamma,omega,phi;/*parameter define*/
    int i,j,k,n1,n2,n3,sn1,pow1,pow2;
    double h;
    float  *inpsg,*xt,*xt1,*re;
    sf_complex *out1;
    sf_file in, out, restor=NULL;
    bool verb,cosine;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");

    if (NULL != sf_getstring ("restor")) {
        restor = sf_input("restor");
    } else {
	restor = NULL;
	sn1=1;
    }
    
   
    sf_setformat(out,"native_complex");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histfloat(in,"d2",&d2)) d2=1;
    if (!sf_histfloat(in,"o2",&o2)) o2=1;
    if (!sf_histint(in,"n3",&n3)) n3=1;
    if (!sf_histfloat(in,"d3",&d3)) d3=1;
    if (!sf_histfloat(in,"o3",&o3)) o3=1;
    
    if(!sf_getfloat("gamma",&gamma)) gamma=0.75;
    /*strength of external force*/
    if(!sf_getfloat("omega",&omega)) omega=1;
    /*angular frequence of external force*/
    if(!sf_getfloat("kxi",&kxi)) kxi=1;
    /*adjustment for input signal*/
    if(!sf_getfloat("x0",&x0)) x0=0;
    /*initial value of x0*/
    if(!sf_getfloat("y0",&y0)) y0=0;
    /*initial value of y0*/
    if(!sf_getint("pow1",&pow1)) pow1=1;
    /*power of first non-linear restitution term*/
    if(!sf_getint("pow2",&pow2)) pow2=3;
    /*power of second non-linear restitution term*/
    if(!sf_getfloat("phi",&phi)) phi=0.;
    /*phase of cosine signal unit=pi*/
    if(!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if(!sf_getbool("cosine",&cosine)) cosine = true;
    /* if n need extenal input for periodic restoring force */
    
    if (NULL != restor){
	
	if (!sf_histint(restor,"n1",&sn1)) sf_error("No n1= in s");
	if (n1>sn1) sf_error("insufficient lenth of external force");
    }
    
    
    h=d1;
    /*step lenth of R-K4*/
    
    re = sf_floatalloc(sn1);
    inpsg = sf_floatalloc(n1*n2);
    xt = sf_floatalloc(n1*n2);
    xt1 = sf_floatalloc(n1*n2);
    out1 = sf_complexalloc(n1*n2);
    
    
    if (NULL != restor) {
	sf_floatread(re,sn1,restor);
    }else {
	for (i=0;i<sn1; i++)
	    re[i] = 0;
    }

    
    if(cosine) {
	for(k = 0; k < n3 ; k++) {
	    sf_floatread(inpsg,n1*n2,in);
	    x00=x0,y00=y0;
	    t0=0;
	    if(verb) sf_warning("step %d of %d;",k,n3);
	    for(j = 0;j < n2;j++) {
		x00=x0,y00=y0;
		t0=0;
		for(i = 0;i < n1;i++){
		    
		    k11 = function1(t0,x00,y00,omega);
		    k12 = function2(t0,x00,y00,inpsg[i+j*n1],kxi,
				    gamma,omega,pow1,pow2,phi);
		    k21 = function1(t0+h/2,x00+k11*h/2,y00+k12*h/2,omega);
		    k22 = function2(t0+h/2,x00+k11*h/2,
				    y00+k12*h/2,inpsg[i+j*n1],kxi,
				    gamma,omega,pow1,pow2,phi);
		    k31 = function1(t0+h/2,x00+k21*h/2,y00+k22*h/2,omega);
		    k32 = function2(t0+h/2,x00+k21*h/2,
				    y00+k22*h/2,inpsg[i+j*n1],kxi,
				    gamma,omega,pow1,pow2,phi);
		    k41 = function1(t0+h,x00+k31*h,y00+k32*h,omega);
		    k42 = function2(t0+h,x00+k31*h,y00+k32*h,inpsg[i+j*n1],
				    kxi,gamma,omega,pow1,pow2,phi);
		    
		    x00+= h*(k11+2*k21+2*k31+k41)/6;
		    y00+= h*(k12+2*k22+2*k32+k42)/6;
		    t0+= h;
		    
		    xt[i+j*n1]=x00,xt1[i+j*n1]=y00;
		}	
	    }
	    for(i=0;i < n1*n2;i++) {
		out1[i] = sf_cmplx(xt[i],xt1[i]);
	    }
	    sf_complexwrite(out1,n1*n2,out);
	}
    }else {
	for(k = 0; k < n3 ; k++) {
	    sf_floatread(inpsg,n1*n2,in);
	    x00=x0,y00=y0;
	    t0=0;
	    if(verb) sf_warning("step %d of %d;",k,n3);
	    for(j = 0;j < n2;j++) {
		x00=x0,y00=y0;
		t0=0;
		for(i = 0;i < n1;i++) {
		    	    
		    k11 = function1(t0,x00,y00,omega);
		    k12 = function3(t0,x00,y00,inpsg[i+j*n1],kxi,
				    gamma,omega,pow1,pow2,re[i]);
		    k21 = function1(t0+h/2,x00+k11*h/2,y00+k12*h/2,omega);
		    k22 = function3(t0+h/2,x00+k11*h/2,
				    y00+k12*h/2,inpsg[i+j*n1],kxi,
				    gamma,omega,pow1,pow2,re[i]);
		    k31 = function1(t0+h/2,x00+k21*h/2,y00+k22*h/2,omega);
		    k32 = function3(t0+h/2,x00+k21*h/2,
				    y00+k22*h/2,inpsg[i+j*n1],kxi,
				    gamma,omega,pow1,pow2,re[i]);
		    k41 = function1(t0+h,x00+k31*h,y00+k32*h,omega);
		    k42 = function3(t0+h,x00+k31*h,y00+k32*h,inpsg[i+j*n1],
				    kxi,gamma,omega,pow1,pow2,re[i]);
		    
		    x00+= h*(k11+2*k21+2*k31+k41)/6;
		    y00+= h*(k12+2*k22+2*k32+k42)/6;
		    t0+= h;
		    
		    
		    xt[i+j*n1]=x00,xt1[i+j*n1]=y00;
		}
	    }
	    for(i=0;i < n1*n2;i++) {
		out1[i] = sf_cmplx(xt[i],xt1[i]);
	    }
	    sf_complexwrite(out1,n1*n2,out);
	}
    }
    
    exit(0);
}

/* 	$Id$	 */
