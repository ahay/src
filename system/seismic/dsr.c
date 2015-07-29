/* DSR modeling/migration in v(z) */
/*
  Copyright (C) 2006 Colorado School of Mines
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

#include "dsr.h"
#include "pshift.h"

static float eps, *vt, *vz, *eta, dw, dz, da;
static int nz, nw, na;
static sf_complex *pp;
static kiss_fftr_cfg forw, invs;

void dsr_init (float eps1 /* regularization */, 
	       int nw1    /* number of frequencies */,
	       int nt     /* time samples */, 
	       float dt   /* time sampling */, 
	       int nz1    /* depth samples */, 
	       float dz1  /* depth sampling */, 
	       float *vt1 /* velocity/slowness */, 
	       float *vz1 /* vertical velocity/slowness */,
	       float *eta1 /* eta */,
	       bool depth /* depth or time migration */,
	       char rule   /* interpolation rule */,
	       int na1     /* angle samples */,
	       float da1   /* angle sampling */)
/*< initialize >*/
{
    eps = eps1;     
    nz = nz1; 
    dz = dz1;
    vt = vt1;
    vz = vz1;
    eta = eta1;
    na = na1;
    da = da1;
    nw = nw1;

    /* determine frequency sampling */
    dw = 2.0*SF_PI/(nt*dt);
    
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
    
    /* allocate workspace */
    pp = sf_complexalloc (nt/2+1);
    pshift_init(depth,rule);
} 

void dsr_close ()
/*< free workspace >*/
{    
    free (pp);  
    free (forw);
    free (invs);
}

static double Power(float a, int p)
{
    int i;
    double b=a;
    for(i=0; i<p-1; ++i)
	b *=a;
    return(b);
}


static float costheta(float v, float vz, float n, float w, float s, float r)
{
    float a;
    double bo2a,coa,a2;

    a2=s;
    s *=s;
    r *=r;

    bo2a=(-4*Power(n,2)*Power(r,2)*Power(s,2) +   
	  2*n*r*s*(r + s)*v*
	  Power(w,2) - r*s*                  
	  (Power(v,2) + Power(vz + 4*n*vz,2))*Power(w,4) +     
	  (1 + 4*n)*(r + s)*v*Power(vz,2)*    
	  Power(w,6) - Power(v,2)*Power(vz,2)*Power(w,8))/     
	((2*n*Power(r,2) + r*(-v + vz + 4*n*vz)*       
	  Power(w,2) - v*vz*Power(w,4))*                     
	 (2*n*Power(s,2) + s*(-v + vz + 4*n*vz)*      
	  Power(w,2) - v*vz*Power(w,4)));
                    
    coa=Power(4*Power(n,2)*Power(r,2)*Power(s,2) -               
	      2*n*r*s*(r + s)*v*
	      Power(w,2) + r*s*                  
	      (Power(v,2) - Power(vz + 4*n*vz,2))*Power(w,4) +     
	      (1 + 4*n)*(r + s)*v*Power(vz,2)*    
	      Power(w,6) - Power(v,2)*Power(vz,2)*Power(w,8),2)/   
	(Power(-2*n*Power(r,2) +                                
	       r*(v - (1 + 4*n)*vz)*Power(w,2) +          
	       v*vz*Power(w,4),2)*                                 
	 Power(-2*n*Power(s,2) +
	       s*(v - (1 + 4*n)*vz)*Power(w,2) +
	       v*vz*Power(w,4),2));
    /*  bo2a=(-1 + (1 + 4*n)*(Power(r,2) + Power(s,2))*Power(vz,2) - 
	2*(1 + 4*n + 8*Power(n,2))*Power(r,2)*Power(s,2)*Power(vz,4) + 
	2*n*Power(r,2)*Power(s,2)*(Power(r,2) + Power(s,2))*Power(vz,6) - 
	4*Power(n,2)*Power(r,4)*Power(s,4)*Power(vz,8))/
	((-1 + 2*n*Power(r,2)*Power(vz,2)*(2 + Power(r,2)*Power(vz,2)))*
	(-1 + 2*n*Power(s,2)*Power(vz,2)*(2 + Power(s,2)*Power(vz,2))));
	coa=Power(1 - (1 + 4*n)*(Power(r,2) + Power(s,2))*Power(vz,2) + 
	8*n*(1 + 2*n)*Power(r,2)*Power(s,2)*Power(vz,4) + 
	2*n*Power(r,2)*Power(s,2)*(Power(r,2) + Power(s,2))*Power(vz,6) - 
	4*Power(n,2)*Power(r,4)*Power(s,4)*Power(vz,8),2)/
	(Power(1 - 2*n*Power(r,2)*Power(vz,2)*(2 + Power(r,2)*Power(vz,2)),2)*
	Power(1 - 2*n*Power(s,2)*Power(vz,2)*(2 + Power(s,2)*Power(vz,2)),2));*/


    if(a2<0)
	a2 = -bo2a -sqrt(SF_ABS(Power(bo2a,2) - coa));
    else
	a2 = -bo2a +sqrt(SF_ABS(Power(bo2a,2) - coa));
    /*warn("bo2a=%f coa=%f a2=%f",bo2a,coa,a2);*/

    if(a2>=0) a=sqrt(a2);
    else a=sqrt(-a2);
  
    a = SF_MIN(1,sqrt(0.5+0.5*a));
    /*kk = sqrt(Power(w,2)*vz-Power(dw*eps,2)*vz-s)+sqrt(Power(w,2)*vz-Power(dw*eps,2)-r);
      w2 = sf_cmplx(dw*eps,w);
      w2 = w2*w2;
      w2 = csqrtf(w2*vz+s)+csqrtf(w2*vz+r);
      a = 0.5/(sqrt(vz)*w)*hypotf(cimagf(w2),kx);*/
    /*warn("a=%f kk=%f w=%f v=%f s2=%f r2=%f",a,kk,w,vz,s,r);*/

    return(a);
}

void dsr (char rule /* rule for angle gathers */,
	  bool inv /* modeling or migration */, 
	  float kx /* midpoint wavenumber */, 
	  float kh /* half-offset wavenumber */, 
	  float *p /* time trace */, 
	  float **q /* angle gathers */)
/*< apply >*/
{
    int iz,iw, ia;
    float s, r, a;
    sf_complex w, k, kk;

    /* convert from midpoint/offset wavenumbers to
       source/receiver wavenumbers */
    s = 0.5*(kx-kh);
    r = 0.5*(kx+kh);
    s *= s;
    r *= r;

    if (inv) { /* modeling */
	for (iw=0; iw<nw; iw++) {
	    pp[iw] = sf_cmplx(q[nz-1][0],0.);
	}
	
	/* loop over migrated times z */
	for (iz=nz-2; iz>=0; iz--) {
	    /* loop over frequencies w */
	    for (iw=0; iw<nw; iw++) {
		w = sf_cmplx(eps*dw,iw*dw);

#ifdef SF_HAS_COMPLEX_H
		k = pshift(w,r,vt[iz],vt[iz+1],vz[iz],eta[iz])+
		    pshift(w,s,vt[iz],vt[iz+1],vz[iz],eta[iz]);
		pp[iw] = q[iz][0] + pp[iw]*cexpf(-k*dz);
#else
		k = sf_cadd(pshift(w,r,vt[iz],vt[iz+1],vz[iz],eta[iz]),
			    pshift(w,s,vt[iz],vt[iz+1],vz[iz],eta[iz]));
		pp[iw] = sf_cadd(sf_cmplx(q[iz][0],0.),
				 sf_cmul(pp[iw],cexpf(sf_crmul(k,-dz))));
#endif
	    }
	}

	kiss_fftri(invs,(const kiss_fft_cpx *) pp, p);
    } else { /* migration */

	/* FFT time -> frequency */
	kiss_fftr(forw, p, (kiss_fft_cpx *) pp);

	/* loop over migrated times z */
	for (iz=0; iz<nz-1; iz++) {
	    /* loop over frequencies w */
	    for (iw=1; iw<nw; iw++) {
		if(iw*dw*iw*dw*vt[iz]/(1+2*eta[iz])<s) continue;
		if(iw*dw*iw*dw*vt[iz]/(1+2*eta[iz])<r) continue;
		w = sf_cmplx(eps*dw,iw*dw);

		/* find angle */

#ifdef SF_HAS_COMPLEX_H
		k = pshift(w,r,vt[iz],vt[iz+1],vz[iz],eta[iz])+
		    pshift(w,s,vt[iz],vt[iz+1],vz[iz],eta[iz]);
		kk = pshift(w,r,vt[iz],vt[iz+1],vt[iz],0.0)+
		    pshift(w,s,vt[iz],vt[iz+1],vt[iz],0.0);
#else
		k = sf_cadd(pshift(w,r,vt[iz],vt[iz+1],vz[iz],eta[iz]),
			    pshift(w,s,vt[iz],vt[iz+1],vz[iz],eta[iz]));
		kk = sf_cadd(pshift(w,r,vt[iz],vt[iz+1],vt[iz],0.0),
			     pshift(w,s,vt[iz],vt[iz+1],vt[iz],0.0));
#endif
		
		if (rule=='a') {
		    a = costheta(vt[iz],vt[iz],eta[iz],cimagf(w),0.5*(kx-kh),0.5*(kx+kh));
		    /* ss = 0.5/(sqrt(vt[iz])*cimagf(w))*hypotf(cimagf(kk),kx); */
		    /*warn("a=%f ss=%f eps=%f w=%f kx=%f s=%f r=%f eta=%f",a,ss,eps,cimagf(w),kx,sqrt(s),sqrt(r),eta[iz]);*/

		    if (a <= 1.) {
			a = acosf(a)/da;
			ia = floorf(a);
			a -= ia;
			if (ia >=0 && ia < na-1) {
			    /* accumulate image (summed over frequency) */
			    q[iz][ia] += (1.-a)*crealf(pp[iw]);
			    q[iz][ia+1] += a*crealf(pp[iw]);
			} else if (ia==na-1) {
			    q[iz][ia] += crealf(pp[iw]);
			}
		    }
		
		} else{
		    a = 0.5/(sqrt(vt[iz])*cimagf(w))*hypotf(cimagf(kk),kx);

		    if (a <= 1.) {
			a = acosf(a)/da;
			ia = floorf(a);
			a -= ia;
			if (ia >=0 && ia < na-1) {
			    /* accumulate image (summed over frequency) */
			    q[iz][ia] += (1.-a)*crealf(pp[iw]);
			    q[iz][ia+1] += a*crealf(pp[iw]);
			} else if (ia==na-1) {
			    q[iz][ia] += crealf(pp[iw]);
			}
		    } 
		}

#ifdef SF_HAS_COMPLEX_H
		pp[iw] *= conjf(cexpf(-k*dz));
#else
		pp[iw] = sf_cmul(pp[iw],conjf(cexpf(sf_crmul(k,-dz))));
#endif
	    }
	}

	for (iw=0; iw<nw; iw++) {
	    for (ia=0; ia < na; ia++) {
		q[nz-1][ia] += crealf(pp[iw]);
	    }
	}
    }
}

/* 	$Id: dsr.c 9238 2012-10-02 02:07:38Z sfomel $	 */

