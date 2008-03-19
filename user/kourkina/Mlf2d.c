/* Convert Dix velocity to interval velocity. 

Input in (x0,t0), output in (x,z).
*/
/*
  Copyright (C) 2008 New York University
  
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <rsf.h>

#include "fastmarch.h"

#define INFTY 1.0e+6

static void init(void);
static void qp(void);

/*************************************************************/

static int nc, neval;
static int nx,nz,nt,nx1,nx2,nz1,nt1,nxt,nxz;
static float hx,hz,ht;
static float *f,*x0,*t0,*v,*q,*s,*g;
static float dx,dz;

/*************************************************************/

/* We will shift the mesh (x,y): [0,xmax-xmin]*[0,ymax-ymin] */

/*****************************************************/

static void init(void) 
{
    int i,j,ind;

    for( i=0;i<nx;i++ ) {
	*(x0+i)=i*hx;
	*(t0+i)=0.0;
	*(v+i)=(*(f+i));	
	*(s+i)=1/(*(f+i));
	*(q+i)=1.0;
	*(g+i)=0.0;
	ind=i;
	for( j=1;j<nz;j++ ) {
	    ind+=nx;
	    *(x0+ind)=0.0;
	    *(t0+ind)=INFTY;
	    *(v+ind)=0.0;
	}
	for( j=1;j<nt;j++ ){
	    ind=i+j*nx;
	    *(q+ind)=1.0;
	    *(g+ind)=0.0;
	}
    }

    fastmarch_init(nx,nz,nt,
		   hx,hz,ht,
		   x0,t0,v,s);
}


/******** COMPUTE Q AND P IN TIME COORDINATES *********/

static void qp(void) 
{
  int i,k=0,ind;
  const float lam=ht/(hx*hx);
  float ff,qq,f0,f1,f2,f3,q0,q1,q2,q3;

  sf_warning("lam=%.4e",lam);

  for( k=0;k < nt1;k++) {
      for( i=1; i<nx1; i++ ) {
	  ind=i+nx*k;
	  q0=*(q+ind-1);
	  q1=*(q+ind+1);
	  q2=(i>1) ? *(q+ind-2) : 1.0;
	  q3=(i<nx2) ? *(q+ind+2) : 1.0;
	  qq=*(q+ind);
	  f0=*(f+ind-1);
	  f1=*(f+ind+1);
	  f2=(i>1) ? *(f+ind-2) : *(f+k*nx);
	  f3=(i<nx1) ? *(f+ind+2) : *(f+nx1+k*nx);
	  ff=*(f+ind);
	  *(g+ind+nx)=0.5*(*(g+ind-1)+(*(g+ind+1)))- 
	      0.25*lam*((f3*q3-ff*qq)/(q1)-(ff*qq-f2*q2)/(q0))/(ff*qq);
	  f2=*(f+ind+nx);
	  *(q+ind+nx)=-1.0/(-1.0/(qq)+0.5*ht*(ff*ff*(*(g+ind))+f2*f2*(*(g+ind+nx))));
      }
  }
  for( i=0; i<nx; i++ ) {
      for( k=0; k<nt; k++ ) {
	  ind=i+nx*k;
	  *(s+ind)=1.0/(*(f+ind)*(*(q+ind)));
      }
  }
}

int main(int argc, char* argv[]) 
{
    bool cheb;
    sf_file fv,fv2,fx,ft;

    sf_init(argc,argv);
    fv = sf_input("in");
    fv2 = sf_output("out");
    fx = sf_output("x0");
    ft = sf_output("t0");
    
    if (SF_FLOAT != sf_gettype(fv)) sf_error("Need float input");
    if (!sf_histint(fv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(fv,"n2",&nt)) sf_error("No n2= in input");

    if (!sf_histfloat(fv,"d1",&hx)) sf_error("No d1= in input");
    if (!sf_histfloat(fv,"d2",&ht)) sf_error("No d2= in input");
    ht /= 2; /* convert to one-way traveltime */

    if (!sf_getint("nz",&nz)) sf_error("Need nz=");
    if (!sf_getfloat("dz",&hz)) sf_error("Need dz=");

    sf_putint(fv2,"n2",nz);
    sf_putfloat(fv2,"d2",hz);
    sf_putfloat(fv2,"o2",0.);

    sf_putint(fx,"n2",nz);
    sf_putfloat(fx,"d2",hz);
    sf_putfloat(fx,"o2",0.);

    sf_putint(ft,"n2",nz);
    sf_putfloat(ft,"d2",hz);
    sf_putfloat(ft,"o2",0.);

    nxt=nx*nt;
    nxz=nx*nz;

    nx1=nx-1;
    nx2=nx-2;

    nz1=nz-1;
    nt1=nt-1;

    dx=0.1*hx;
    dz=0.1*hz;
 
    f= sf_floatalloc(nxt);
    x0= sf_floatalloc(nxz);
    t0= sf_floatalloc(nxz);
    v= sf_floatalloc(nxz);
    s= sf_floatalloc(nxt); 
    q= sf_floatalloc(nxt);
    g= sf_floatalloc(nxt);

    sf_floatread(f,nxt,fv);
    
    /* three parts of the algorithm */
    init();
    qp();
    fastmarch();

    sf_floatwrite(v,nxz,fv2);
    sf_floatwrite(x0,nxz,fx);
    sf_floatwrite(t0,nxz,ft);
    
    exit(0);
}  
	   			
