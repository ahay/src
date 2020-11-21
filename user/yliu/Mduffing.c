/* Duffing differential equation solved by 4th order Runge-Kutta method. 
Duffing equation: x'' + 0.5 x' - x + x^3 = gamma cos(omega t) + kxi input(t)
*/
/*
  Copyright (C) 2013 Jilin University
  
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

float function2(float t,float x,float y,float Rd,float kxi,float gamma,float omega,int pow1,int pow2)
{
  float y1;
  y1=omega*(powf(x,pow1)-powf(x,pow2)-0.5*y+kxi*Rd+gamma*cos(omega*t));
  return y1;
}

float function3(float t,float x,float y,float Rd,float kxi,float gamma,float omega,int pow1,int pow2,float s)
{
  float y1;
  y1=omega*(powf(x,pow1)-powf(x,pow2)-0.5*y+kxi*Rd+gamma*s);
  return y1;
}

int main(int argc,char *argv[])
{
    float x0,y0,t0,d1,o1,k11,k21,k31,k41,k12,k22,k32,k42;
    float kxi,gamma,omega;/*parameter define*/
    int i,n1,n2,sn1,pow1,pow2,i2;
    double h;
    float  *rd,*xt,*xt1,*s;
    sf_complex *out1;
    sf_file in, out, sfile=NULL;
    bool verb,ricker;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");

    if (NULL != sf_getstring ("sfile"))
    {
        sfile = sf_input("sfile");
    } else {
	sfile = NULL;
    }
    

    sf_setformat(out,"native_complex");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    /*check n1*/
    n2 = sf_leftsize(in,1);

    if (NULL != sfile){
	if (!sf_histint(sfile,"n1",&sn1)) sf_error("No n1= in s");
	if (n1>sn1) sf_error("insufficient lenth of external force");
    }

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
    if(!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if(!sf_getbool("ricker",&ricker)) ricker = false;
    /* if y need extenal input for external force */
    
    h=d1;
    /*step lenth of R-K4*/
    
    s = sf_floatalloc(n1);
    rd = sf_floatalloc(n1);
    xt = sf_floatalloc(n1);
    xt1 = sf_floatalloc(n1);
    out1 = sf_complexalloc(n1);

    for(i2=0; i2 < n2; i2++) {

	sf_floatread(rd,n1,in);
	
	if (NULL != sfile){
	    sf_floatread(s,n1,sfile);
	}else{
	    for (i=0;i<n1; i++)
		s[i] = 0;
	}
	
	t0=0;
	if(ricker){
	    for(i = 0;i < n1;i++)
	    {
		if(verb) sf_warning("step %d t0=%f;",i,t0);
		
		k11 = function1(t0,x0,y0,omega);
		k12 = function3(t0,x0,y0,rd[i],kxi,gamma,omega,pow1,pow2,s[i]);
		k21 = function1(t0+h/2,x0+k11*h/2,y0+k12*h/2,omega);
		k22 = function3(t0+h/2,x0+k11*h/2,y0+k12*h/2,rd[i],
				kxi,gamma,omega,pow1,pow2,s[i]);
		k31 = function1(t0+h/2,x0+k21*h/2,y0+k22*h/2,omega);
		k32 = function3(t0+h/2,x0+k21*h/2,y0+k22*h/2,rd[i],
				kxi,gamma,omega,pow1,pow2,s[i]);
		k41 = function1(t0+h,x0+k31*h,y0+k32*h,omega);
		k42 = function3(t0+h,x0+k31*h,y0+k32*h,rd[i],kxi,
				gamma,omega,pow1,pow2,s[i]);
		
		x0+= h*(k11+2*k21+2*k31+k41)/6;
		y0+= h*(k12+2*k22+2*k32+k42)/6;
		t0+= h;
		;
		
		xt[i]=x0,xt1[i]=y0;
	    }
	}else{
	    for(i = 0;i < n1;i++)
	    {
		
		if(verb) sf_warning("step %d t0=%f;",i,t0);
		
		k11 = function1(t0,x0,y0,omega);
		k12 = function2(t0,x0,y0,rd[i],kxi,gamma,omega,pow1,pow2);
		k21 = function1(t0+h/2,x0+k11*h/2,y0+k12*h/2,omega);
		k22 = function2(t0+h/2,x0+k11*h/2,y0+k12*h/2,rd[i],
				kxi,gamma,omega,pow1,pow2);
		k31 = function1(t0+h/2,x0+k21*h/2,y0+k22*h/2,omega);
		k32 = function2(t0+h/2,x0+k21*h/2,y0+k22*h/2,rd[i],
				kxi,gamma,omega,pow1,pow2);
		k41 = function1(t0+h,x0+k31*h,y0+k32*h,omega);
		k42 = function2(t0+h,x0+k31*h,y0+k32*h,rd[i],kxi,
				gamma,omega,pow1,pow2);
		
		x0+= h*(k11+2*k21+2*k31+k41)/6;
		y0+= h*(k12+2*k22+2*k32+k42)/6;
		t0+= h;
		
		xt[i]=x0,xt1[i]=y0;
	    }
	}
	
	
	for(i=0;i < n1;i++){
	    out1[i] = sf_cmplx(xt[i],xt1[i]);
	}
	sf_complexwrite(out1,n1,out);
    }

    exit(0);
}

/* 	$Id$	 */
