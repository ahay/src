/* 2D/3D Velocity analysis by using Duffing differential equation solved by 4th order Runge-Kutta method. 
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
    int i,j,k,m,n1,sn1,n3,p,a,b,lx,pow1,pow2,*w;
    int iw,vn,t0n,xn,midp,winsz,starp;
    float x0,y0,x00,y00,dut0,h,k11,k21,k31,k41,k12,k22,k32,k42;
    float kxi,gamma,omega,delta,phi;
    float dt,dx,t0,t00,o1,deltat0,v0,v00,dv,f,gx;
    float *orig,*output,*winds,*xt,*yt,*re;
    char time[]="Time",vel[]="Velocity",unit2[]="s",unit1[]="m/s";
    bool verb,cosine;
    sf_file cmp,outf,restor=NULL;;
    
    sf_init (argc, argv); 
    cmp = sf_input("in");
    outf = sf_output("out");
    
    if (NULL != sf_getstring ("restor")) {
        restor = sf_input("restor");
    } else {
	restor = NULL;
	sn1=1;
    }
    
    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histint(cmp,"n2",&xn)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");
    n3 = sf_leftsize(cmp,2);
    
    if(!sf_getint("winsz",&winsz)) winsz=200;
    /* for each trace,the width of window.unit:ms*/
    if(!sf_getfloat("v0",&v0)) v0=1000;
    /* init Vel for velocity scan */
    if(!sf_getfloat("dv",&dv)) dv=20;
    /* step lenth for velocity scan */
    if(!sf_getint("vn",&vn)) vn=100;
    /* numbers of velscan*/
    if(!sf_getfloat("t0",&t0)) t0=o1;
    /* t0 scan start point */
    if(!sf_getfloat("deltat0",&deltat0)) deltat0=dt;
    /* step lenth for t0 scan */
    if(!sf_getint("t0n",&t0n)) t0n=n1;
    /* numbers of t0scan*/
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
    if(!sf_getbool("cosine",&cosine)) cosine = true;
    /* if n need extenal input for periodic restoring force */
    if (!sf_getfloat("delta",&delta)) delta=0.01;
    /*The density of judgement grid*/
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if(!sf_getfloat("gx",&gx)) gx=2.0;
    /*Size of grid*/
    
    if (NULL != restor){
	if (!sf_histint(restor,"n1",&sn1)) sf_error("No n1= in s");
	if (winsz*xn>sn1) sf_error("insufficient lenth of external force");
    }
    
    lx = (fabs(gx)*2)/delta;
    
    orig = sf_floatalloc(n1*xn);
    winds = sf_floatalloc(winsz*xn);
    re = sf_floatalloc(sn1);
    xt = sf_floatalloc(winsz*xn);
    yt = sf_floatalloc(winsz*xn);
    output = sf_floatalloc(vn*t0n);
    w = sf_intalloc(lx*lx);
    
    sf_putint(outf, "n1", vn);
    sf_putfloat(outf, "d1", dv);
    sf_putfloat(outf, "o1", v0);
    sf_putint(outf, "n2", t0n);
    sf_putfloat(outf, "d2", deltat0);
    sf_putfloat(outf, "o2", t0);
    sf_putstring(outf,"label1",vel);
    sf_putstring(outf,"label2",time);
    sf_putstring(outf,"unit1",unit1);
    sf_putstring(outf,"unit2",unit2);
    for(m=0; m < n3; m++) {
	if(verb) sf_warning("n3step = %d of %d;\n",m,n3-1);
	sf_floatread(orig,n1*xn,cmp);
	
	h=dt;
	/*step lenth of R-K4*/
	t00=t0;
	for(k = 0; k < t0n; k++) {
	    v00 = v0;
	    if(verb) sf_warning("step = %d of %d;",k,t0n);
	    for(j = 0; j < vn; j++) {
		
		for(i = 0; i < xn; i++){
		    f = t00*t00+powf(i*dx,2)/powf(v00,2);
		    midp = sqrtf(f)/dt;
		    starp = midp-(winsz/2);
		    
		    for(iw = 0; iw < winsz; iw++) {
			if(((starp+iw) > n1)||(starp < 0)) {
			    winds[iw + winsz*i]= 0 ;
			}else{
			    winds[iw + winsz*i]=orig[starp+iw + i*n1];
			}
		    }		
		}
		if(cosine){
		    x00=x0,y00=y0;
		    dut0=0;
		    for(i = 0;i < winsz*xn;i++){
			k11 = function1(dut0,x00,y00,omega);
			k12 = function2(dut0,x00,y00,winds[i],kxi,
					gamma,omega,pow1,pow2,phi);
			k21 = function1(dut0+h/2,x00+k11*h/2,
					y00+k12*h/2,omega);
			k22 = function2(dut0+h/2,x00+k11*h/2,
					y00+k12*h/2,winds[i],kxi,
					gamma,omega,pow1,pow2,phi);
			k31 = function1(dut0+h/2,x00+k21*h/2,
					y00+k22*h/2,omega);
			k32 = function2(dut0+h/2,x00+k21*h/2,
					y00+k22*h/2,winds[i],kxi,
					gamma,omega,pow1,pow2,phi);
			k41 = function1(dut0+h,x00+k31*h,y00+k32*h,omega);
			k42 = function2(dut0+h,x00+k31*h,y00+k32*h,winds[i],
					kxi,gamma,omega,pow1,pow2,phi);
			
			x00+= h*(k11+2*k21+2*k31+k41)/6;
			y00+= h*(k12+2*k22+2*k32+k42)/6;
			dut0+= h;
			
			xt[i]=x00,yt[i]=y00;
		    }
		} else {
		    x00=x0,y00=y0;
		    dut0=0;
		    for(i = 0;i < winsz*xn;i++) {	    
			k11 = function1(dut0,x00,y00,omega);
			k12 = function3(dut0,x00,y00,winds[i],kxi,
					gamma,omega,pow1,pow2,re[i]);
			k21 = function1(dut0+h/2,x00+k11*h/2,
					y00+k12*h/2,omega);
			k22 = function3(dut0+h/2,x00+k11*h/2,
					y00+k12*h/2,winds[i],kxi,
					gamma,omega,pow1,pow2,re[i]);
			k31 = function1(dut0+h/2,x00+k21*h/2,
					y00+k22*h/2,omega);
			k32 = function3(dut0+h/2,x00+k21*h/2,
					y00+k22*h/2,winds[i],kxi,
					gamma,omega,pow1,pow2,re[i]);
			k41 = function1(dut0+h,x00+k31*h,y00+k32*h,omega);
			k42 = function3(dut0+h,x00+k31*h,y00+k32*h,winds[i],
					kxi,gamma,omega,pow1,pow2,re[i]);
			
			x00+= h*(k11+2*k21+2*k31+k41)/6;
			y00+= h*(k12+2*k22+2*k32+k42)/6;
			dut0+= h;
			
			xt[i]=x00,yt[i]=y00;
			
		    }
		}
		
		for(i=0;i < lx*lx; i++){
		    w[i]=0;
		}
		
		for(i=0;i<winsz*xn;i++){
		    a = (xt[i]+gx)/delta;
		    b = (yt[i]+gx)/delta;
		    w[b*lx+a]=1;
		}
		
		p = 0;
		for(i = 0; i<lx*lx; i++){
		    p=p+w[i];
		}
		output[j+k*vn] = p;
		
		v00 += dv;
	    }
	    t00+= deltat0;
	    
	}
	sf_floatwrite(output,vn*t0n,outf);
    }
    exit(0);
}

/* 	$Id: Mduffing2.c 11371 2013-11-21 03:03:21Z yang_liu $	 */
