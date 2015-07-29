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
#include "chebyshev.h"

#define INFTY 1.0e+6
#define NC 1000

static void init(char meth);
static void chebyshev_init(void);
static void qp_lf(void);
static void qp_cheb(void);

/*************************************************************/

static int nc, neval;
static int nx,nz,nt,nx1,nx2,nz1,nt1,nxt,nxz;
static float hx,hz,ht;
static float *f,*x0,*t0,*v,*q,*s,*g,*fc,*y;
static float dx,dz,xmax;
static float cp[NC];
static float bma,bpa;
static const float eps=0.;

/*************************************************************/

/* We will shift the mesh (x,y): [0,xmax-xmin]*[0,ymax-ymin] */

/*****************************************************/

static void init(char meth) 
{
    int i,j,ind;

    for( i=0;i<nx;i++ ) {
	*(x0+i)=i*hx;
	*(t0+i)=0.0;
	*(v+i)=(*(f+i));	
	*(s+i)=1/(*(f+i));
	*(q+i)=1.0;
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
	}
    }

    if ('c'==meth) chebyshev_init();
    fastmarch_init(nx,nz,nt,
		   hx,hz,ht,
		   x0,t0,v,s);
}

static void chebyshev_init(void)
{
    int i,j,ind,k,klo,khi;
    float yp1,ypn, *u, *b, *b2,a1,a2,un,/*qn,*/x,aux;

    u = sf_floatalloc(nx);
    b = sf_floatalloc(nx);
    b2 = sf_floatalloc(nx);

    /*--- compute f in chebyshev points---*/
    for( i=nc-1; i>=0; i-- ) {
	cp[nc-i-1]=cosf(SF_PI*(i+0.5)/nc);
    }
    for(j=0;j<nt;j++) {
	ind=j*nx;
	/* compute cubicspline coefficients */
	yp1=(*(f+1+ind)-(*(f+ind)))/hx;
	ypn=(*(f+ind+nx1)-(*(f+ind+nx2)))/hx;
	for(i=0;i<nx;i++) b[i]=*(f+ind+i);
	b2[0]=-0.5;
	u[0]=(3.0/hx)*((b[1]-b[0])/hx-yp1);
	for( i=1; i<nx-1; i++ ) {
	    aux=0.5*b2[i-1]+2.0;
	    b2[i]=-0.5/aux;
	    u[i]=(b[i+1]-b[i])/hx-(b[i]-b[i-1])/hx;
	    u[i]=(3.0*u[i]/hx-0.5*u[i-1])/aux;
	}
/*	qn=0.5; */
	un=(3.0/hx)*(ypn-(b[nx1]-b[nx2])/hx);
	b2[nx1]=(un-0.5*u[nx2])/(0.5*b2[nx2]+1.0);
	for( k=nx-2; k>=0; k-- ) b2[k]=b2[k]*b2[k+1]+u[k];
	/* evaluate f at chebyshev points */
	for(i=0;i<nc;i++) {
	    x=cp[i]*bma+bpa;
	    klo=0;
	    khi=nx1;
	    while(khi-klo>1) {
		k=(khi+klo)>>1;
		if(k*hx>x) khi=k;
		else klo=k;
	    }
	    a1=(khi*hx-x)/hx;
	    a2=(x-klo*hx)/hx;
	    *(fc+j*nc+i)=a1*b[klo]+a2*b[khi]+
		((a1*a1*a1-a1)*b2[klo]+(a2*a2*a2-a2)*b2[khi])*hx*hx/6.0;
	}
    }

    free(u);
    free(b);
    free(b2);
}

/******** COMPUTE Q AND P IN TIME COORDINATES *********/

static void qp_lf(void) 
/* Lax-Friedrichs algorithm */
{
  int i,k=0,ind;
  const float lam=ht/(hx*hx);
  float ff,qq,/* f0,f1, */ f2,f3,q0,q1,q2,q3;

  sf_warning("lam=%.4e",lam);

  for( i=0;i<nxt;i++ ) {
      *(g+i)=0.0;
  }

  for( k=0;k < nt1;k++) {
      for( i=1; i<nx1; i++ ) {
	  ind=i+nx*k;
	  q0=*(q+ind-1);
	  q1=*(q+ind+1);
	  q2=(i>1) ? *(q+ind-2) : 1.0;
	  q3=(i<nx2) ? *(q+ind+2) : 1.0;
	  qq=*(q+ind);
/*	  f0=*(f+ind-1);
	  f1=*(f+ind+1); */
	  f2=(i>1) ? *(f+ind-2) : *(f+k*nx);
	  f3=(i<nx1) ? *(f+ind+2) : *(f+nx1+k*nx);
	  ff=*(f+ind);
	  *(g+ind+nx)=0.5*(*(g+ind-1)+(*(g+ind+1)))- 
	      0.25*lam*((f3*q3-ff*qq)/(q1)-(ff*qq-f2*q2)/(q0))/(ff*qq);
	  f2=*(f+ind+nx);
	  *(q+ind+nx)=-1.0/
	      (-1.0/(qq)+0.5*ht*(ff*ff*(*(g+ind))+f2*f2*(*(g+ind+nx))));
      }
  }
  for( i=0; i<nx; i++ ) {
      for( k=0; k<nt; k++ ) {
	  ind=i+nx*k;
	  *(s+ind)=1.0/(*(f+ind)*(*(q+ind)));
      }
  }
}

static void qp_cheb(void) 
/* Chebyshev spectral algorithm */
{
  int i,k,ind,ind1;
  float b[NC],b2[NC],coef[NC],cder[NC],rhs[3][NC],creg[NC],crder[NC],crd2[NC];
  float f0,f1,yp1,ypn;
  const float a1=23.0/12.0, a2=-4.0/3.0, a3=5.0/12.0;
  char ch='y';

  for(i=0;i<nc;i++) {
    ind=i;
    *(y+ind)=-1.0;
    *(g+ind)=0.0;
  }
  for( k=0;k<2;k++ ) {
    /*---chebspectrum for y---*/
    for(i=0;i<nc;i++) b[i]=*(y+i+k*nc);
    chebcoef(b,coef,nc);
    /*---computation of RHS---*/
    for(i=0;i<nc;i++) {
      ind=i+k*nc;
      b[i]=(*(fc+ind))/(*(y+ind));
    }
    chebcoef(b,coef,nc);
    for(i=0;i<nc;i++) cder[i]=0.0;
    chebder(0.,xmax,coef,cder,nc);
    for(i=0;i<nc;i++) {
      ind=i+nc*k;
      b[i]=chebeval(cder,cp[i],neval)*(*(y+ind));
    }
    chebcoef(b,coef,nc);
    chebder(0.,xmax,coef,cder,nc);
    /*---add (y_t/f^2)_{xx}---*/
    for(i=0;i<nc;i++) {
      ind=i+nc*k;
      b[i]=(*(g+ind));
    }
    chebcoef(b,creg,nc);
    for(i=0;i<nc;i++) {crder[i]=0.0; crd2[i]=0.0;}
    chebder(0.,xmax,creg,crder,nc);
    chebder(0.,xmax,crder,crd2,nc);
    /*---time step---*/
    for(i=0;i<nc;i++) {
      ind=i+k*nc;
      b[i]=chebeval(cder,cp[i],neval);
      rhs[k][i]=(*(y+ind))*b[i]/(*(fc+ind))+eps*chebeval(crd2,cp[i],neval);
      ind1=ind+nc;
      *(g+ind1)=*(g+ind)+ht*rhs[k][i];
      f0=*(fc+ind);
      f1=*(fc+ind1);
      *(y+ind1)=*(y+ind)+0.5*ht*(f0*f0*(*(g+ind))+f1*f1*(*(g+ind1)));
    }
  }
  while(k<nt1  && ch=='y') {
      sf_warning("%d of %d;",k,nt1);

    /*---chebspectrum for y---*/
    for(i=0;i<nc;i++) b[i]=*(y+i+k*nc);
    chebcoef(b,coef,nc);
    /*---computation of RHS---*/
    for(i=0;i<nc;i++) {
      ind=i+k*nc;
      b[i]=(*(fc+ind))/(*(y+ind));
    }
    chebcoef(b,coef,nc);
    chebder(0.,xmax,coef,cder,nc);
    for(i=0;i<nc;i++) {
      ind=i+nc*k;
      b[i]=chebeval(cder,cp[i],neval)*(*(y+ind));
    }
    chebcoef(b,coef,nc);
    chebder(0.,xmax,coef,cder,nc);
    for(i=0;i<nc;i++) {
      ind=i+nc*k;
      b[i]=(*(g+ind));
    }
    chebcoef(b,creg,nc);
    for(i=0;i<nc;i++) {crder[i]=0.0; crd2[i]=0.0;}
    chebder(0.,xmax,creg,crder,nc);
    chebder(0.,xmax,crder,crd2,nc);
    /*---time step---*/
    for(i=0;i<nc;i++) {
      ind=i+k*nc;
      b[i]=chebeval(cder,cp[i],neval);
      rhs[k%3][i]=(*(y+ind))*b[i]/(*(fc+ind))+eps*chebeval(crd2,cp[i],neval);
      ind1=ind+nc;
      *(g+ind1)=*(g+ind)+
	  ht*(a1*rhs[k%3][i]+a2*rhs[(k-1)%3][i]+a3*rhs[(k-2)%3][i]);
      f0=*(fc+ind);
      f1=*(fc+ind1);
      *(y+ind1)=*(y+ind)+0.5*ht*(f0*f0*(*(g+ind))+f1*f1*(*(g+ind1)));
      if( *(y+ind1)>=-0.1 ) ch='n';
    }
    k++;
  }
  sf_warning(".");
  /*---------- find the velocity on the regular mesh by cubicsplines --------*/
  for( k=1;k<nt;k++) {
      for(i=0;i<nc;i++) b[i]=*(y+i+nc*k);
      yp1=(*(y+k*nc+1)-(*(y+k*nc)))/((cp[1]-cp[0]));
      ypn=(*(y+k*nc+nc-1)-(*(y+k*nc+nc-2)))/((cp[nc-1]-cp[nc-2]));
      spline(cp,b,nc,yp1,ypn,b2);
      for( i=0;i<nx;i++ ) {
	  ind=i+nx*k;
	  q[ind]=-1.0/splineval(cp,b,b2,nc,((i*hx)-bpa)/bma);
	  s[ind]=1.0/(f[ind]*q[ind]);
      }
  }
}

int main(int argc, char* argv[]) 
{
    float xmin;
    char *method, meth;
    sf_file fv=NULL,fv2=NULL,fx=NULL,ft=NULL;

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

    if (NULL == (method = sf_getstring("method"))) method="lf";
    /* method (chebyshev,lax-friedrichs) */
    meth = method[0];

    if ('c'==meth) {
	if (!sf_getint("nc",&nc)) nc=100;
	/* number of chebyshev coefficients */
	if (nc > NC) sf_error("nc must be smaller than %d",NC);
	if (!sf_getint("neval",&neval)) neval=20;
	/* numvber of used chebyshev coefficients */
    }

    nxt=nx*nt;
    nxz=nx*nz;

    nx1=nx-1;
    nx2=nx-2;

    nz1=nz-1;
    nt1=nt-1;

    xmax = nx1*hx;
    xmin = 0.;

    bma=0.5*(xmax-xmin);
    bpa=0.5*(xmax+xmin);

    dx=0.1*hx;
    dz=0.1*hz;
 
    f= sf_floatalloc(nxt);
    x0= sf_floatalloc(nxz);
    t0= sf_floatalloc(nxz);
    v= sf_floatalloc(nxz);
    s= sf_floatalloc(nxt); 
    q= sf_floatalloc(nxt);

    if ('c'==meth) {
	g= sf_floatalloc(nt*NC);
	fc= sf_floatalloc(nt*NC);
	y= sf_floatalloc(nt*NC);
    } else {
	g= sf_floatalloc(nxt);
	fc= NULL;
	y= NULL;
    } 

    sf_floatread(f,nxt,fv);

    /* three parts of the algorithm */
    init(meth);

    if ('c'==meth) {
	qp_cheb();
    } else {
	qp_lf();
    }

    fastmarch();

    sf_floatwrite(v,nxz,fv2);
    sf_floatwrite(x0,nxz,fx);
    sf_floatwrite(t0,nxz,ft);


    exit(0);
}
