/* VTI eikonal solver. */
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
/*^*/

#include "wkbjti.h"

#define LHD 20
#define NHD 1+2*LHD

#define TINY 1.0e-4
#define CFL 0.98
#define EPSS 1.0e-8

void eikpex (int na, float da, float r, float dr, 
	float sc[], float uc[], float wc[], float tc[],
	     float sn[], float un[], float wn[], float tn[]);
void sigma (int na, float da, float r, float dr, 
	float uc[], float wc[], float sc[],
	float un[], float wn[], float sn[]);
void beta (int na, float da, float r, float dr, 
	float uc[], float wc[], float bc[],
	float un[], float wn[], float bn[]);
void tripp(int n, float *d, float *e, float *c, float *b);
float svp(float a, float a1111, float a3333,float a1133,float a1313);
float dsvp(float a, float a1111, float a3333,float a1133,float a1313);
float dsvg(float a, float a1111, float a3333,float a1133,float a1313, float *angp);
void recttopolar (
	int nx, float dx, float fx, int ny, float dy, float fy, float **p,
	int na, float da, float fa, int nr, float dr, float fr, float **q);
void polartorect (
	int na, float da, float fa, int nr, float dr, float fr, float **q,
	int nx, float dx, float fx, int ny, float dy, float fy, float **p);


/* Numerically solving the differential source linear equation */
void wkbjti (float xs                /* x source location*/, 
		float zs                  /* z source location */, 
	        int nz                   /* z dimension */, 
	        float dz                   /*  z sampling interval*/, 
		float fz                    /* z first sample */, 
	        int nx                   /* x dimension */, 
	        float dx                   /*  x sampling interval*/, 
		float fx                    /* x first sample */, 
		float** a1111               /*  a1111 */,
		float** a3333              /* a3333 */,
		float** a1313              /* a1313  */,
	        float** a1133              /* a1133  */,
		float** time            /* time */,
	     float** angle            /* angle */,
	     float** sig            /* sig */,
	     float** bet            /* bet */)
/*< Run VTI eikonal solver >*/
{
	int ia,ir,na,nr;
	float a,r,da,dr,fa,fr,ex,ez,ea,rmax,rmaxs,uu,wor,angp,
		**sp,**tp,**up,**wp,**ap;
	float **a1111p,**a3333p,**a1313p,**a1133p;
	

	/* shift coordinates so source is at (x=0,z=0) */
	fx -= xs;
	fz -= zs;
	ex = fx+(nx-1)*dx;
	ez = fz+(nz-1)*dz;
	/* determine polar coordinate sampling */
	rmaxs = fx*fx+fz*fz;
	rmaxs = SF_MAX(rmaxs,fx*fx+ez*ez);
	rmaxs = SF_MAX(rmaxs,ex*ex+ez*ez);
	rmaxs = SF_MAX(rmaxs,ex*ex+fz*fz);
	rmax = sqrt(rmaxs);
	dr = SF_MIN(SF_ABS(dx),SF_ABS(dz));
	nr = 1+SF_NINT(rmax/dr);
	dr = rmax/(nr-1);
	fr = 0.0;
	if (fx==0.0 && fz==0.0) {
		fa = 0.0;  ea = SF_PI/2.0;
	} else if (fx<0.0 && fz==0.0) {
		fa = -SF_PI/2.0;  ea = SF_PI/2.0;
	} else if (fx==0.0 && fz<0.0) {
		fa = 0.0;  ea = SF_PI;
	} else {
		fa = -SF_PI;  ea = SF_PI;
	}
	da = dr/rmax;
	na = 1+SF_NINT((ea-fa)/da);
	da = (ea-fa)/(na-1);
	if (fa==-SF_PI && ea==SF_PI)
		na = na-1;
	/* allocate space */
	a1111p = sf_floatalloc2(na,nr);
	a3333p = sf_floatalloc2(na,nr);
	a1313p = sf_floatalloc2(na,nr);
	a1133p = sf_floatalloc2(na,nr);
	sp = sf_floatalloc2(na,nr);
	tp = sf_floatalloc2(na,nr);
	up = sf_floatalloc2(na,nr);
	wp = sf_floatalloc2(na,nr);
	ap = sf_floatalloc2(na,nr);
	/*sf_error("na=%d da=%f nr=%d dr=%f nx=%d dx=%f dz=%f nz=%d vh=%f ea=%f fa=%f",na,da,nr,dr,nx,dx,dz,nz,sqrt(a1111[0][0]),ea,fa);*/
	/* convert from rectangular to polar coordinates */
	recttopolar(nz,dz,fz,nx,dx,fx,a1111,na,da,fa,nr,dr,fr,a1111p);
	recttopolar(nz,dz,fz,nx,dx,fx,a3333,na,da,fa,nr,dr,fr,a3333p);
	recttopolar(nz,dz,fz,nx,dx,fx,a1133,na,da,fa,nr,dr,fr,a1133p);
	recttopolar(nz,dz,fz,nx,dx,fx,a1313,na,da,fa,nr,dr,fr,a1313p);

	/* compute traveltimes and derivatives in source region */
	for (ia=0, a=fa; ia<na; ++ia,a+=da) {
	  tp[0][ia] = 0.0;
	  up[0][ia] = sp[0][ia] = svp(a,a1111p[0][ia],a3333p[0][ia],
				      a1133p[0][ia],a1313p[0][ia]);
	  wp[0][ia] = 0.0;
	  ap[0][ia] = a;
	}

	for (ia=0, a=fa; ia<na; ++ia,a+=da)
	  /*tp[1][ia] = dr*svg(a,a1111p[0][ia],a3333p[0][ia],
			     a1133p[0][ia],a1313p[0][ia]);*/
	  /*tp[1][ia] = dr*sqrt(a3333p[0][ia]);*/
	  tp[1][ia] = dr*svp(a,.5*(a1111p[0][ia]+a1111p[1][ia]),.5*(a3333p[0][ia]+a3333p[1][ia]),
			     .5*(a1133p[0][ia]+a1133p[1][ia]),.5*(a1313p[0][ia]+a1313p[1][ia]));


	for (ia=0, a=fa; ia<na; ++ia,a+=da) {
	  /*wp[1][ia] = -.5*(tp[1][ia+1]-tp[1][ia-1])/da;*/
	  /*t11 = dr*svp(a-.01,.5*(a1111p[0][ia]+a1111p[1][ia]),.5*(a3333p[0][ia]+a3333p[1][ia]),
			     .5*(a1133p[0][ia]+a1133p[1][ia]),.5*(a1313p[0][ia]+a1313p[1][ia]));
	  t12 = dr*svp(a+.01,.5*(a1111p[0][ia]+a1111p[1][ia]),.5*(a3333p[0][ia]+a3333p[1][ia]),
			     .5*(a1133p[0][ia]+a1133p[1][ia]),.5*(a1313p[0][ia]+a1313p[1][ia]));*/
	  /*t11 = dr*svp(a-.01,.5*(a1111p[0][ia]),.5*(a3333p[0][ia]),
			     .5*(a1133p[0][ia]),.5*(a1313p[0][ia]));
	  t12 = dr*svp(a+.01,.5*(a1111p[0][ia]),.5*(a3333p[0][ia]),
			     .5*(a1133p[0][ia]),.5*(a1313p[0][ia]));*/
	  wp[1][ia] = SF_SIG(a)*dr*dsvg(a,.5*(a1111p[0][ia]+a1111p[1][ia]),.5*(a3333p[0][ia]+a3333p[1][ia]),
			     .5*(a1133p[0][ia]+a1133p[1][ia]),.5*(a1313p[0][ia]+a1313p[1][ia]),&angp);
	  /*wp[1][ia] = SF_SIG(a)*dr*dsvg(a,a1111p[1][ia],a3333p[1][ia],
				      a1133p[1][ia],a1313p[1][ia],&angp);*/
	  /*wp[1][ia] = 0;*/
	  wp[1][ia] = SF_SIG(a)*dr*dsvp(a,.5*(a1111p[0][ia]+a1133p[1][ia]),.5*(a3333p[0][ia]+a3333p[1][ia]),
			     .5*(a1133p[0][ia]+a1133p[1][ia]),.5*(a1313p[0][ia]+a1313p[1][ia]));
	  sp[1][ia] = svp(angp,a1111p[1][ia],a3333p[1][ia],
				      a1133p[1][ia],a1313p[1][ia]);
	  /*sp[1][ia] = 1/sqrt(a3333p[1][ia]);*/
	  wor=wp[1][ia]/dr;
	  uu = sp[1][ia]*sp[1][ia]-wor*wor;
	  if(uu<=0) sf_error("\tRaypath has a too large curvature!\n\t A smoother velocity is required. \n");
	  up[1][ia] = sqrt(uu);
	  /*fprintf(stderr,"ia=%d a=%f angp=%f wp=%f sp=%f uu=%f up=%f\n",ia,a,angp,wp[1][ia],sp[1][ia],uu,up[1][ia]);*/
	  ap[1][ia] = a+asin(wp[1][ia]/(sp[1][ia]*dr));
	  /*fprintf(stderr,"ap=%f angp=%f\n",ap[1][ia],angp);*/

	}		
	/*wp[1][0] = 2*wp[1][1]-wp[1][2];
	wp[1][na-1] = 2*wp[1][na-2]-wp[1][na-3];
	sp[1][0] = svp(fa,a1111p[1][0],a3333p[1][0],
				      a1133p[1][0],a1313p[1][0]);
	sp[1][na-1] = svp(fa+(na-1)*da,a1111p[1][na-1],a3333p[1][na-1],
				      a1133p[1][na-1],a1313p[1][na-1]);
	wor   = wp[1][0]/dr;
	uu    = sp[1][0]*sp[1][0]-wor*wor;
	if(uu<=0) err("\tRaypath has a too large curvature!\n\t A smoother velocity is required. \n");
	up[1][0] = sqrt(uu);
	wor   = wp[1][na-1]/dr;
	uu    = sp[1][na-1]*sp[1][na-1]-wor*wor;
	if(uu<=0) err("\tRaypath has a too large curvature!\n\t A smoother velocity is required. \n");
	up[1][na-1] = sqrt(uu);
	ap[1][0] = fa+asin(wp[1][0]/(sp[1][0]*dr));
	ap[1][na-1] = fa+(na-1)*da+asin(wp[1][na-1]/(sp[1][na-1]*dr));*/

/*	for (ia=0; ia<na; ++ia,a+=da)
	  fprintf(stderr,"tp=%f wp=%f sp=%f up=%f ap=%f\n",tp[1][ia],wp[1][ia],sp[1][ia],up[1][ia],ap[1][ia]);*/

/* 	tt=cpusec();   */

	/* solve eikonal equation for remaining times and derivatives */
	for (ir=1,r=dr; ir<nr-1; ++ir,r+=dr) {
	        for (ia=0; ia<na; ++ia)
		    sp[ir+1][ia] = svp(ap[ir][ia],a1111p[ir+1][ia],a3333p[ir+1][ia],
				       a1133p[ir+1][ia],a1313p[ir+1][ia]);
		eikpex(na,da,r,dr,
			sp[ir],up[ir],wp[ir],tp[ir],
			sp[ir+1],up[ir+1],wp[ir+1],tp[ir+1]);
		/*fprintf(stderr,"sp=%f tp=%f ap=%f\n",sp[ir][200],tp[ir][200],ap[ir][200]);*/
		for (ia=0,a=fa; ia<na; ++ia,a+=da)
		   ap[ir+1][ia] = a+asin(wp[ir+1][ia]/(sp[ir+1][ia]*(r+dr)));
	}
	
	/* convert times from polar to rectangular coordinates */
	polartorect(na,da,fa,nr,dr,fr,tp,nz,dz,fz,nx,dx,fx,time);

/*  	fprintf(stderr,"\t CPU time for traveltimes= %f \n",cpusec()-tt); 
 	tt=cpusec();   */
	
	polartorect(na,da,fa,nr,dr,fr,ap,nz,dz,fz,nx,dx,fx,angle);
/*  	fprintf(stderr,"\t CPU time for propagation angles= %f\n", 	
		cpusec()-tt); 
	tt=cpusec();   */
	
	/* compute sigmas  for initial values */
	for (ir=0,r=0; ir<2; ++ir,r+=dr) 
		for (ia=0; ia<na; ++ia) tp[ir][ia] = r/sp[ir][ia];

	/* solve diffrence equation for remaining sigmas */
	for (ir=1,r=dr; ir<nr-1; ++ir,r+=dr) 
 		sigma(na,da,r,dr,up[ir],wp[ir],tp[ir],
			up[ir+1],wp[ir+1],tp[ir+1]);  
	polartorect(na,da,fa,nr,dr,fr,tp,nz,dz,fz,nx,dx,fx,sig);

/* 	fprintf(stderr,"\t CPU time for sigmas= %f \n",cpusec()-tt); 
	tt=cpusec(); */
	
	/* compute betas for initial values */
	for (ir=0; ir<2; ++ir) 
		for (ia=0,a=fa; ia<na; ++ia,a+=da) tp[ir][ia] = a;

	/* solve diffrence equation for remaining betas */
	for (ir=1,r=dr; ir<nr-1; ++ir,r+=dr) 
 		beta(na,da,r,dr,up[ir],wp[ir],tp[ir],
			up[ir+1],wp[ir+1],tp[ir+1]);  
	polartorect(na,da,fa,nr,dr,fr,tp,nz,dz,fz,nx,dx,fx,bet);
	
/* 	fprintf(stderr,"\t CPU time for incident angles= %f \n",
		cpusec()-tt); */

	free(a1111p);
	free(a3333p);
	free(a1133p);
	free(a1313p);
	free(sp);
	free(tp);
	free(up);
	free(wp);
	free(ap);
	

}


void eikpex (int na, float da, float r, float dr, 
	float sc[], float uc[], float wc[], float tc[],
	float sn[], float un[], float wn[], float tn[])
/*****************************************************************************
Eikonal equation extrapolation of times and derivatives in polar coordinates
******************************************************************************
Input:
na		number of a samples
da		a sampling interval
r		current radial distance r
dr		radial distance to extrapolate
sc		array[na] of slownesses at current r
uc		array[na] of dt/dr at current r
wc		array[na] of dt/da at current r
tc		array[na] of times t at current r
sn		array[na] of slownesses at next r

Output:
un		array[na] of dt/dr at next r (may be equivalenced to uc)
wn		array[na] of dt/da at next r (may be equivalenced to wc)
tn		array[na] of times t at next r (may be equivalenced to tc)
******************************************************************************
Notes:
If na*da==2*PI, then the angular coordinate is wrapped around (periodic). 

This function implements the finite-difference method described by Bill
Symes (Rice University) and Jos van Trier (Stanford University) in a
(1990) preprint of a paper submitted to Geophysics.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 07/16/90
******************************************************************************/
{
	int i,wrap;
	float drleft,drorig,frac,cmax,umaxl,uminr,uminm,umaxm,oda=1/da,
		uu,unew,uold,ueol,ueor,wor,invr,*wtemp,*s;
	float a,b,c,x;
	/*float xx[na],duda[na][4];
	
	for (i=0; i<na; ++i)
	  xx[i] =i*da;*/
	/* allocate workspace */
	wtemp = sf_floatalloc(na);
	s = sf_floatalloc(na);
	
	/* remember the step size */
	drleft = drorig = dr;
	
	/* initialize slownesses to values at current r */
	for (i=0; i<na; ++i)
		s[i] = sc[i];
	
	/* copy inputs to output */
	for (i=0; i<na; ++i) {
		un[i] = uc[i];
		wn[i] = wc[i];
		tn[i] = tc[i];
	}
	
	/* determine if angular coordinate wraps around */
	wrap = SF_ABS(na*da-2.0*SF_PI)<0.01*SF_ABS(da);
	
	/* loop over intermediate steps with adaptive stepsize */
	while (drleft>0.0) {
		
		/* determine adaptive step size according to CFL condition */
		for (i=0,cmax=TINY; i<na; ++i) {
			if (r*SF_ABS(un[i])<TINY*SF_ABS(wn[i]))
				cmax = 1.0/TINY;
			else
				cmax = SF_MAX(cmax,SF_ABS(wn[i]/(r*un[i])));
		}
		dr = SF_MIN(drleft,CFL/cmax*r*da);

		if(dr<drleft)
		  fprintf(stderr,"r=%f dr=%f drleft=%f\n",r,dr,drleft);
		
		/* if angles wrap around */
		if (wrap) {
			umaxl = (wn[na-1]>0.0 ? un[na-1] : s[0]);
			if (wn[0]>0.0) {
				uminm = s[0];
				umaxm = un[0];
			} else {
				uminm = un[0];
				umaxm = s[0];
			}
			uminr = (wn[1]>0.0 ? s[0] : un[1]);
			ueol = uminm+umaxl;
			ueor = uminr+umaxm;
			wtemp[0] = wn[0]+dr*(ueor-ueol)*oda;
			umaxl = (wn[na-2]>0.0 ? un[na-2] : s[na-1]);
			if (wn[na-1]>0.0) {
				uminm = s[na-1];
				umaxm = un[na-1];
			} else {
				uminm = un[na-1];
				umaxm = s[na-1];
			}
			uminr = (wn[0]>0.0 ? s[na-1] : un[0]);
			ueol = uminm+umaxl;
			ueor = uminr+umaxm;
			wtemp[na-1] = wn[na-1]+dr*(ueor-ueol)*oda;
		
		/* else, if angles do not wrap around */
		} else {
			if (wn[0]<=0.0)
				wtemp[0] = wn[0] + 
					dr*(un[1]-un[0])*oda; 
			else
				wtemp[0] = 0.0;
			if (wn[na-1]>=0.0) 
				wtemp[na-1] = wn[na-1] +
					dr*(un[na-1]-un[na-2])*oda;
			else
				wtemp[na-1] = 0.0;
		}
		
		/* update interior w values via Enquist/Osher scheme */
		for (i=1; i<na-1; ++i) {
		        /*if((wn[i]<0 && wn[i-1]>0 && wn[i+1]<0) || 
			   (wn[i]>0 && wn[i-1]<0 && wn[i+1]>0)) {
			  fprintf(stderr,"r=%f ia-1=%d wn=%f s=%f un=%f\n",r,i-1,wn[i-1],s[i-1],un[i-1]);
			  fprintf(stderr,"r=%f ia=%d wn=%f s=%f un=%f\n",r,i,wn[i],s[i],un[i]);
			  fprintf(stderr,"r=%f ia+1=%d wn=%f s=%f un=%f\n",r,i,wn[i+1],s[i+1],un[i+1]);
			}*/
			umaxl = (wn[i-1]>0.0 ? un[i-1] : s[i]);
			if (wn[i]>0.0) {
				uminm = s[i];
				umaxm = un[i];
			} else {
				uminm = un[i];
				umaxm = s[i];
			}
			uminr = (wn[i+1]>0.0 ? s[i] : un[i+1]);
			ueol = uminm+umaxl;
			ueor = uminr+umaxm;
			wtemp[i] = wn[i]+dr*(ueor-ueol)*oda;
		}

		/*for (i=1; i<na-1; ++i) {
		  if(ABS(wn[i])<0.00*r*s[i])
		    wtemp[i] = wn[i]+.5*dr*(un[i+1]-un[i-1])*oda;
		  else {
		    if(wn[i]>0)
		      wtemp[i] = wn[i]+dr*(un[i]-un[i-1])*oda;
		    else
		      wtemp[i] = wn[i]+dr*(un[i+1]-un[i])*oda;
		  }
		}*/

		/*fprintf(stderr,"wtemp=%f\n",wtemp[160]);*/
		/*daa= dr/drorig*da;
		for (i=1; i<na-1; ++i) {
		  a = 0.5*(un[i+1]-2*un[i]+un[i-1])*oda*oda;
		  b = 0.5*(un[i+1]-un[i-1])*oda;
		  if(wn[i]>0)
		    wtemp[i] = wn[i]+dr*(-a*daa+b);
		  else
		    wtemp[i] = wn[i]+dr*(a*daa+b);
		}
		fprintf(stderr,"wtemp=%f\n",wtemp[160]);*/
		
		/*for (i=1; i<na-1; ++i) {
		  a = 0.5*(un[i+1]-2*un[i]+un[i-1])*oda*oda;
		  b = 0.5*(un[i+1]-un[i-1])*oda;
		  if(wn[i]>0)
		    wtemp[i] = wn[i]+dr*(-a*daa+b);
		  else
		    wtemp[i] = wn[i]+dr*(a*daa+b);
		}*/

		for (i=1; i<na-1; ++i) {
		  /*b = SF_SIG(wn[i])*pow(ABS(wn[i]/(r*s[i])),.3);
		  b = SF_SIG(wn[i]);*/
		  x = wn[i]/(r*s[i]);
		  b = 2*x/(1+x*x);
		  b = SF_SIG(x)*(1-pow(SF_ABS(x)-1,200));
		  /*b = x;
		  b = SF_SIG(x);*/
		  a = 0.5*(1-b);
		  c = 1-a; 
		  wtemp[i] = wn[i]+dr*(a*un[i+1]+b*un[i]-c*un[i-1])*oda;
		}
/*fprintf(stderr,"wtemp=%f b=%f c=%f a=%f x=%f wn=%f s=%f\n",wtemp[160],b,c,a,x,wn[160],s[160]);*/

		/*cakima (na,xx,un,duda);
		for (i=0; i<na; ++i)
		  wtemp[i] = wn[i]+dr*duda[i][1];*/
		
		
		/* decrement the size of step left to do */
		drleft -= dr;
		
		/* update radial coordinate and its inverse */
		r += dr;
		invr = 1.0/r;
		
		/* linearly interpolate slowness for new r */
		frac = drleft/drorig;
		for (i=0; i<na; ++i)
			s[i] = frac*sc[i]+(1.0-frac)*sn[i];
		
		/* update w and u; integrate u to get t */
		for (i=0; i<na; i++) {
			wn[i] = wtemp[i];
			wor = wn[i]*invr;
			uu = (s[i]-wor)*(s[i]+wor);
			/*if(uu<=0) err("\tRaypath has a too large curvature!\n\t A smoother velocity is required. \n");*/
			uu=MAX(0,uu);
 			unew = sqrt(uu); 
			uold = un[i];
			un[i] = unew;
			tn[i] += 0.5*dr*(unew+uold);
		}
	}
	
	/* free workspace */
	free(wtemp);
	free(s);
}

void sigma (int na, float da, float r, float dr, 
	float uc[], float wc[], float sc[],
	float un[], float wn[], float sn[])
/*****************************************************************************
difference equation extrapolation of "sigma" in polar coordinates
******************************************************************************
Input:
na		number of a samples
da		a sampling interval
r		current radial distance r
dr		radial distance to extrapolate
uc		array[na] of dt/dr at current r
wc		array[na] of dt/da at current r
sc		array[na] of sigma  at current r
un		array[na] of dt/dr at next r
wn		array[na] of dt/da at next r

Output:
sn		array[na] of sigma at next r 
******************************************************************************

This function implements the Crank-Nicolson finite-difference method with
boundary conditions dsigma/da=0.
******************************************************************************
Author:  Zhenyue Liu, Colorado School of Mines, 07/8/92
******************************************************************************/
{
	int i;
	float r1,*d,*b,*c,*e;
	
	/* allocate workspace */
	d = sf_floatalloc(na-2);
	b = sf_floatalloc(na-2);
	c = sf_floatalloc(na-2);
	e = sf_floatalloc(na-2);
	
	r1 = r+dr;
 	
	/* Crank-Nicolson */
 	for (i=0; i<na-2; ++i) {
		d[i] = (uc[i+1]+un[i+1])/(2.0*dr);
		e[i] = (wn[i+1]/(r1*r1)+wc[i+1]/(r*r))/(8.0*da);
		b[i] = 1.0-(sc[i+2]-sc[i])*e[i]
			+d[i]*sc[i+1];
		c[i] = -e[i];
	} 
	d[0] += c[0];
	d[na-3] += e[na-3]; 
	
	tripp(na-2,d,e,c,b);
	for(i=0;i<na-2; ++i) sn[i+1]=b[i];
	sn[0] = sn[1];
	sn[na-1] = sn[na-2];
	
	
	/* free workspace */
	free(d);
	free(c);
	free(e);
	free(b);
}

void beta (int na, float da, float r, float dr, 
	float uc[], float wc[], float bc[],
	float un[], float wn[], float bn[])
/*****************************************************************************
difference equation extrapolation of "beta" in polar coordinates
******************************************************************************
Input:
na		number of a samples
da		a sampling interval
r		current radial distance r
dr		radial distance to extrapolate
uc		array[na] of dt/dr at current r
wc		array[na] of dt/da at current r
bc		array[na] of beta  at current r
un		array[na] of dt/dr at next r
wn		array[na] of dt/da at next r

Output:
bn		array[na] of beta at next r 
******************************************************************************

This function implements the Crank-Nicolson finite-difference method, with 
boundary conditions dbeta/da=1. 
******************************************************************************
author:  Zhenyue Liu, Colorado School of Mines, 07/8/92
******************************************************************************/
{
	int i;
	float r1,*d,*b,*c,*e;
	
	/* allocate workspace */
	d = sf_floatalloc(na-2);
	b = sf_floatalloc(na-2);
	c = sf_floatalloc(na-2);
	e = sf_floatalloc(na-2);
	
	r1 = r+dr;
	/* Crank-Nicolson */
   	for (i=0; i<na-2; ++i) {
		d[i] = uc[i+1]*r*r+un[i+1]*r1*r1;
		e[i] = (wn[i+1]+wc[i+1])*dr/(4.0*da);
		b[i] = -(bc[i+2]-bc[i])*e[i]
			+d[i]*bc[i+1];
		c[i] = -e[i];
	}   
	d[0] += c[0];
	d[na-3] += e[na-3]; 
	b[0] += da*c[0];
	b[na-3] -= da*e[na-3];
	
	tripp(na-2,d,e,c,b);
	for(i=0;i<na-2; ++i) bn[i+1]=b[i];
	bn[0] = bn[1]-da;
	bn[na-1] = bn[na-2]+da;
	
	
	/* free workspace */
	free(d);
	free(c);
	free(e);
	free(b);
}

static void exch(float x, float y);
static void exch(float x, float y)
{    
	float t;
	t=x; x=y; y=t;
}
void tripp(int n, float *d, float *e, float *c, float *b)
/*******************************************************************
Solve an unsymmetric tridiagonal system that uses Gaussian elimination 
with partial pivoting
********************************************************************
Input:
d	diagonal vector of matrix
e       upper-diagonal vector of matrix
c       lower-diagonal vector of matrix
b       right-hand vector
n       dimension of matrix

Output:
b       solution vector
*******************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 7/06/92
*********************************************************************/
{
	int k;
	float temp;

	
/*      elimination   */
	for(k=0; k<n-1; ++k){
	    c[k] = 0;
 	    if(SF_ABS(d[k])<SF_ABS(c[k+1])){
	        exch(d[k],c[k+1]);
		exch(e[k],d[k+1]);
		exch(c[k],e[k+1]);
		exch(b[k],b[k+1]);
		} 
		
	    if(d[k]==0 ) sf_error("coefficient matrix is singular!\n");
	    temp = c[k+1]/d[k];  
	    d[k+1] -= temp*e[k];
	    e[k+1] -= temp*c[k];
	    b[k+1] -= temp*b[k];
        } 
	 
/*      substitution      */
	if(d[n-1]==0 ) sf_error("coefficient matrix is singular!\n");
	b[n-1] = b[n-1]/d[n-1];
	b[n-2] = (b[n-2] - b[n-1]*e[n-2])/d[n-2];		
	for(k=n-3; k>=0; --k)
	    b[k] = (b[k] - b[k+1]*e[k] - b[k+2]*c[k])/d[k];
	    
}	



float svp(float a, float a1111, float a3333,float a1133,float a1313)
{
	float sint,cost,a2,b2,sin2t,cos2t,psi,gamma,sqgamma,v,e,d;

	a2= a3333;
	b2= a1313;
	e = .5*(a1111-a3333)/a3333;
	d = .5*((a1133+a1313)*(a1133+a1313)-(a3333-a1313)*(a3333-a1313))/(a3333*(a3333-a1313));
	sint=sin(a);
	cost=cos(a);
	sin2t=sint*sint;
	cos2t=cost*cost;
	psi=1.-b2/a2;
	gamma=0.25*psi*psi+(2*d-e)*psi*sin2t*cos2t+(psi+e)*e*sin2t*sin2t;
	sqgamma=sqrt(MAX(gamma,0));
	v=sqrt(MAX(a3333*(1.+e*sin2t-0.5*psi+sqgamma),0));

	return(1/v);
}

float dsvp(float a, float a1111, float a3333,float a1133,float a1313)
{
	float sint,cost,a2,b2,sin2t,cos2t,psi,gamma,sqgamma,v,e,d,dgamma,dv;

	a2= a3333;
	b2= a1313;
	e = .5*(a1111-a3333)/a3333;
	d = .5*((a1133+a1313)*(a1133+a1313)-(a3333-a1313)*(a3333-a1313))/(a3333*(a3333-a1313));
	sint=sin(a);
	cost=cos(a);
	sin2t=sint*sint;
	cos2t=cost*cost;
	psi=1.-b2/a2;
	gamma=0.25*psi*psi+(2*d-e)*psi*sin2t*cos2t+(psi+e)*e*sin2t*sin2t;
	dgamma=(2*d-e)*psi*2*sint*cost*cos2t +
             ( 4*(psi+e)*e - 2*(2*d-e)*psi )*sint*sin2t*cost;
	sqgamma=sqrt(gamma);
	v=sqrt(a3333*(1.+e*sin2t-0.5*psi+sqgamma));
	dv=0.5*a2*(2*e*sint*cost + 0.5*dgamma/sqgamma)/v;

	return(-dv/(v*v));
}

float svg(float a, float a1111, float a3333,float a1133,float a1313)
{

  int itr=0,ntr=LHD;
  float sint,cost,a2,b2,sin2t,cos2t,psi,gamma,sqgamma,v,e,d;
  float dv,dgamma,sina,f1,f2,sint1,sint2,vg;
  float tol=TINY,dsin,err,bb,sinda,cosda;

  a = SF_ABS(a); 
  a2= a3333;
  b2= a1313;
  e = .5*(a1111-a3333)/a3333;
  d = .5*((a1133+a1313)*(a1133+a1313)-(a3333-a1313)*(a3333-a1313))/(a3333*(a3333-a1313));
  sint= sina = sin(a);
  cost=sqrt(1-sint*sint);
  sin2t=sint*sint;
  cos2t=cost*cost;
  psi=1.-b2/a2;
  gamma=0.25*psi*psi+(2*d-e)*psi*sin2t*cos2t+(psi+e)*e*sin2t*sin2t;
  dgamma=(2*d-e)*psi*2*sint*cost*cos2t +
             ( 4*(psi+e)*e - 2*(2*d-e)*psi )*sint*sin2t*cost;
  sqgamma=sqrt(gamma);
  v=sqrt(a3333*(1.+e*sin2t-0.5*psi+sqgamma));
  dv=0.5*a2*(2*e*sint*cost + 0.5*dgamma/sqgamma)/v;
  bb=dv/v;
  cosda=1/sqrt(1+bb*bb);
  sinda=bb*cosda;
  f1=sina-sint*cosda-cost*sinda;
  err=SF_ABS(f1);
  if(sint<.5)
    sint1=sint+.01;
  else
    sint1=sint-.01;

  while(err>tol && itr<ntr) {
     cost=sqrt(1-sint1*sint1);
     sin2t=sint1*sint1;
     cos2t=cost*cost;
     psi=1.-b2/a2;
     gamma=0.25*psi*psi+(2*d-e)*psi*sin2t*cos2t+(psi+e)*e*sin2t*sin2t;
     dgamma=(2*d-e)*psi*2*sint1*cost*cos2t +
             ( 4*(psi+e)*e - 2*(2*d-e)*psi )*sint1*sin2t*cost;
     sqgamma=sqrt(gamma);
     v=sqrt(a3333*(1.+e*sin2t-0.5*psi+sqgamma));
     dv=0.5*a2*(2*e*sint1*cost + 0.5*dgamma/sqgamma)/v;
     bb=dv/v;
     cosda=1/sqrt(1+bb*bb);
     sinda=bb*cosda;
     f2=sina-sint1*cosda-cost*sinda;
     err=SF_ABS(f2);
     sint2=sint1;
     dsin=-f2*(sint1-sint)/(f2-f1);
     sint1+=dsin;
     if(sint1>1.0) sint1=1.0;
     if(sint1<0.0) sint1=0.0;
     sint=sint2;
     f1  = f2;
     ++itr;
  }

  if(itr==ntr)
    fprintf(stderr,"no convergence in svg with dsin=%f, sin=%f and err=%f\n",
	    dsin,sint,err);
  
  vg=sqrt(v*v+dv*dv);
  return(1/vg);
}

float dsvg(float a, float a1111, float a3333,float a1133,float a1313, float *angp)
{

  int itr=0,ntr=LHD;
  double sint,cost,a2,b2,sin2t,cos2t,psi,gamma,sqgamma,v,e,d;
  double dv,dgamma,sina,f1,f2,sint1,sint2,vg,dvg;
  double tol=EPSS,dsin,err,bb,sinda,cosda;

  a = SF_ABS(a); 
  a2= a3333;
  b2= a1313;
  e = .5*(a1111-a3333)/a3333;
  d = .5*((a1133+a1313)*(a1133+a1313)-(a3333-a1313)*(a3333-a1313))/(a3333*(a3333-a1313));
  sint= sina = sin(a);
  cost=sqrt(1-sint*sint);
  sin2t=sint*sint;
  cos2t=cost*cost;
  psi=1.-b2/a2;
  gamma=0.25*psi*psi+(2*d-e)*psi*sin2t*cos2t+(psi+e)*e*sin2t*sin2t;
  dgamma=(2*d-e)*psi*2*sint*cost*cos2t +
             ( 4*(psi+e)*e - 2*(2*d-e)*psi )*sint*sin2t*cost;
  sqgamma=sqrt(gamma);
  v=sqrt(a3333*(1.+e*sin2t-0.5*psi+sqgamma));
  dv=0.5*a2*(2*e*sint*cost + 0.5*dgamma/sqgamma)/v;
  bb=dv/v;
  cosda=1/sqrt(1+bb*bb);
  sinda=bb*cosda;
  f1=sina-sint*cosda-cost*sinda;
  err=SF_ABS(f1);
  if(SF_ABS(sint)<.5)
    sint1=sint+SF_SIG(sint)*.01;
  else
    sint1=sint-SF_SIG(sint)*.01;

  while(err>tol && itr<ntr) {
     cost=sqrt(1-sint1*sint1);
     sin2t=sint1*sint1;
     cos2t=cost*cost;
     psi=1.-b2/a2;
     gamma=0.25*psi*psi+(2*d-e)*psi*sin2t*cos2t+(psi+e)*e*sin2t*sin2t;
     dgamma=(2*d-e)*psi*2*sint1*cost*cos2t +
             ( 4*(psi+e)*e - 2*(2*d-e)*psi )*sint1*sin2t*cost;
     sqgamma=sqrt(gamma);
     v=sqrt(a3333*(1.+e*sin2t-0.5*psi+sqgamma));
     dv=0.5*a2*(2*e*sint1*cost + 0.5*dgamma/sqgamma)/v;
     bb=dv/v;
     cosda=1/sqrt(1+bb*bb);
     sinda=bb*cosda;
     f2=sina-sint1*cosda-cost*sinda;
     err=SF_ABS(f2);
     sint2=sint1;
     dsin=-f2*(sint1-sint)/(f2-f1);
     sint1+=dsin;
     if(sint1>1.0) sint1=1.0;
     if(sint1<0.0) sint1=0.0;
     sint=sint2;
     f1  = f2;
     ++itr;
  }

  if(itr==ntr)
    fprintf(stderr,"no convergence in svg with dsin=%f, sin=%f and err=%f\n",
	    dsin,sint,err);
  
  *angp= asin(sint);
  vg=sqrt(v*v+dv*dv);
  /*ddgamma=(2*d-e)*psi*2*(cos2t*cos2t-3*sin2t*cos2t)+
           ( 4*(psi+e)*e - 2*(2*d-e)*psi )*(3*cos2t*sin2t-sin2t*sin2t);
  ddv= 0.5*a2*(2*e*(cos2t-sin2t)+0.5*ddgamma/sqgamma-0.25*dgamma*dgamma/(sqgamma*gamma))/v-
            dv/(v*v);
  dvg=(dv*ddv+v*dv)/vg;*/
  dvg=vg*dv/v;
  return(-dvg/(vg*vg));
}

void recttopolar (
	int nx, float dx, float fx, int ny, float dy, float fy, float **p,
	int na, float da, float fa, int nr, float dr, float fr, float **q)
/*****************************************************************************
Convert a function of p(x,y) to q(a,r), where x = r*cos(a) and y = r*sin(a)
******************************************************************************
Input:
nx		number of x samples
dx		x sampling interval
fx		first x sample
ny		number of y samples
dy		y sampling interval
fy		first y sample
p		array[ny][nx] containing samples of p(x,y)
na		number of a samples
da		a sampling interval
fa		first a sample
nr		number of r samples
dr		r sampling interval
fr		first r sample

Output:
q		array[nr][na] containing samples of q(a,r)
******************************************************************************
Notes:
The polar angle a is measured in radians.

Linear extrapolation is used to determine the value of p(x,y) for
x and y coordinates not in the range corresponding to nx, dx, ....
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/15/90
******************************************************************************/
{
	int ia,ir,ix,iy;
	float a,r,x,y,xi,yi,sx,sy;
	
	/* for all r */
	for (ir=0,r=fr; ir<nr; ++ir,r+=dr) {
	
		/* for all a */
		for (ia=0,a=fa; ia<na; ++ia,a+=da) {
		
			/* determine x and y */
			x = r*cos(a);
			y = r*sin(a);
			
			/* determine sample indices */
			xi = (x-fx)/dx;
			ix = xi;
			if (ix<0 ) xi = ix = 0; 
			if (ix>nx-2) {ix = nx-2; xi = nx-1;}
			yi = (y-fy)/dy;
			iy = yi;
			if (iy<0) yi = iy = 0;
			if (iy>ny-2) {iy = ny-2; yi = ny-1;}
			
			/* bilinear interpolation */
			sx = xi-ix;
			sy = yi-iy;
			q[ir][ia] = (1.0-sy)*((1.0-sx)*p[iy][ix] + 
						sx*p[iy][ix+1]) +
					sy*((1.0-sx)*p[iy+1][ix] +
						sx*p[iy+1][ix+1]);
		}
	}
}

void polartorect (
	int na, float da, float fa, int nr, float dr, float fr, float **q,
	int nx, float dx, float fx, int ny, float dy, float fy, float **p)
/*****************************************************************************
Convert a function of q(a,r) to p(x,y), where x = r*cos(a) and y = r*sin(a)
******************************************************************************
Input:
na		number of a samples
da		a sampling interval
fa		first a sample
nr		number of r samples
dr		r sampling interval
fr		first r sample
nx		number of x samples
dx		x sampling interval
fx		first x sample
ny		number of y samples
dy		y sampling interval
fy		first y sample
q		array[nr][na] containing samples of q(a,r)

Output:
p		array[ny][nx] containing samples of p(x,y)
******************************************************************************
Notes:
The polar angle a is measured in radians.

Linear extrapolation is used to determine the value of q(a,r) for
a and r coordinates not in the range corresponding to na, da, ....
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/15/90
******************************************************************************/
{
	int ix,iy,ia,ir;
	float x,y,a=0.0,r,ai,ri,sa,sr;
	
	/* for all y */
	for (iy=0,y=fy; iy<ny; ++iy,y+=dy) {
	
		/* for all x */
		for (ix=0,x=fx; ix<nx; ++ix,x+=dx) {
		
			/* determine a and r */
			if (x !=0.0)
				a = atan2((double) y,(double) x);
			else if (y>0.0)
				a = SF_PI/2.0;
			else if (y<0.0)
				a = -SF_PI/2.0;
			else if (y==0.0)
				a = 0.0;
			
			r = sqrt(x*x+y*y);
			
			/* determine sample indices */
			ai = (a-fa)/da;
			ia = ai;
			if (ia<0) ai = ia = 0;
			if (ia>na-2) {ai = na-1; ia = na-2;}
			ri = (r-fr)/dr;
			ir = ri;
			if (ir<0) ri = ir = 0;
			if (ir>nr-2) {ri = nr-1; ir = nr-2;}
			
			/* bilinear interpolation */
			sa = ai-ia;
			sr = ri-ir;
			p[iy][ix] = (1.0-sr)*((1.0-sa)*q[ir][ia] + 
						sa*q[ir][ia+1]) +
					sr*((1.0-sa)*q[ir+1][ia] +
						sa*q[ir+1][ia+1]);
		}
	}
}
