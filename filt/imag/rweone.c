/* one-way Rienmannian Wavefield Extrapolation */
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

#include "fft1.h"

#include "rweone.h"

static axa ax,az,aw,ar;
static int method;

static float ***m0, ***n0, ***r0;
static float   *mu,   *nu,   *ro;

static complex float *lk,*rk,*ck;
static complex float *a,*b,*c,*u;

static float c1,c2,sixth;
static float kmu,knu,kro;

static complex float *dax;

static float *tap;
static int   ntap;

static kiss_fft_cfg forw, invs;
static float ffts;
static float okx,dkx;

complex float *wfl,*wtt;
float *mtt,*msk;
int nloop;

#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);

void rweone_init(
    axa ax_,
    axa az_,
    axa aw_,
    axa ar_,
    int method_)
/*< initialize >*/
{
    ax = ax_;
    az = az_;
    aw = aw_;
    ar = ar_;

    method = method_;
    
    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;
    
    /* space for F-D */
    m0 = sf_floatalloc3(az.n,ax.n,ar.n);
    n0 = sf_floatalloc3(az.n,ax.n,ar.n);
    r0 = sf_floatalloc3(az.n,ax.n,ar.n);
    
    mu = sf_floatalloc(ax.n);
    nu = sf_floatalloc(ax.n);
    ro = sf_floatalloc(ax.n);

    if(! sf_getfloat("c1",&c1)) c1=0.47;
    if(! sf_getfloat("c2",&c2)) c2=0.37;
    if(! sf_getfloat("sixth",&sixth)) sixth=0.122996;

    lk = sf_complexalloc(ax.n);
    rk = sf_complexalloc(ax.n);
    ck = sf_complexalloc(ax.n);

    a = sf_complexalloc(ax.n);
    b = sf_complexalloc(ax.n-1);
    c = sf_complexalloc(ax.n-1);
    u = sf_complexalloc(ax.n);

    kmu = 1./ az.d;             /* 1/ dz */
    knu = 1 /(   2*ax.d*ax.d);  /* 1/( 2 dx^2) */
    kro = 1./(az.d*ax.d*ax.d);  /* 1/(dz dx^2) */

    /* data slice */
    dax = sf_complexalloc(ax.n);

    /* taper */
    if(! sf_getint("ntap",&ntap)) ntap=1;
    tap = sf_floatalloc(ax.n);
    rweone_tap_init();

    /* Fourier-domain */
    rweone_fft_init(ax.n);
    okx = -  SF_PI /        ax.d ;
    dkx = 2.*SF_PI /(ax.n * ax.d);

    /* MRS temp arrays */
    wfl = sf_complexalloc(ax.n);
    wtt = sf_complexalloc(ax.n);

    msk = sf_floatalloc(ax.n);
    mtt = sf_floatalloc(ax.n);

    if(! sf_getint("nloop",&nloop)) nloop=3;
    rweone_mrs_init();
}

void rweone_main(
    complex float **dat,
    float         **img,
    float         ** aa,
    float         ** bb,
    float         ** mm,
    float         ** a0,
    float         ** b0)
/*< run one-way extrapolation >*/
{
    int iw,iz,ix,ir;
    float w;

    switch(method) {
	case 3: /* PSC */
	    sf_warning("compute PSC coef");
	    rweone_psc_coef(aa,bb,a0,b0);
	    break;
	case 2: /* FFD */
	    sf_warning("compute FFD coef");
	    rweone_ffd_coef(aa,bb,a0,b0);
	    break;
	case 1: /* SSF */
	    break;
	case 0: /* XFD */
	    sf_warning("compute XFD coef");
	    rweone_xfd_coef(aa,bb);
	    break;
    }

    switch(method) {
	case 3:
	case 2:
	case 1:
	    for(iw=0;iw<aw.n;iw++) {
		w=aw.o+iw*aw.d; w*=2;
		if(method==3)      { sf_warning("PSC %d %d",iw,aw.n); }
		else if(method==2) { sf_warning("FFD %d %d",iw,aw.n); }
		else               { sf_warning("SSF %d %d",iw,aw.n); }
		
		for(ix=0;ix<ax.n;ix++) {
		    dax[ix] = dat[iw][ix];
		}
		
		for(iz=0;iz<az.n;iz++) {
		    for(ix=0;ix<ax.n;ix++) {
			wtt[ix] = dax[ix];
			dax[ix] = 0.;
		    }
		    
		    for(ir=0;ir<ar.n;ir++) {
			for(ix=0;ix<ax.n;ix++) {
			    wfl[ix] = wtt[ix];
			}
			
			rweone_phs(wfl,w,iz,ir,a0,b0);
			rweone_ssf(wfl,w,iz,ir,aa,a0);
			if(method!=1) {
			    rweone_fds(wfl,w,iz,ir);
			}
			rweone_mrs(wfl,  iz,ir,mm,dax);
		    }
		    rweone_tap(dax);
		    for(ix=0;ix<ax.n;ix++) {
			img[ix][iz] += crealf(dax[ix]);
		    }
		}
	    }
	    break;
	case 0:
	default:
	    for(iw=0;iw<aw.n;iw++) {
		w=aw.o+iw*aw.d; w*=2;
		sf_warning("XFD %d %d",iw,aw.n);
		
		for(ix=0;ix<ax.n;ix++) {
		    dax[ix] = dat[iw][ix];
		}
		
		for(iz=0;iz<az.n;iz++) {
		    rweone_ssh(dax,w,iz,aa);
		    rweone_fds(dax,w,iz, 0);
		    rweone_tap(dax);
		    
		    for(ix=0;ix<ax.n;ix++) {
			img[ix][iz] += crealf(dax[ix]);
		    }
		}
	    }
	    break;
    } /* end switch method */
}

/*------------------------------------------------------------*/
/* 
 * Compute coefficients
 */
/*------------------------------------------------------------*/

void rweone_xfd_coef(
    float **aa,
    float **bb)
/*< XFD coefficients >*/
{
    int iz,ix;

    for(ix=0;ix<ax.n;ix++) {
	for(iz=0;iz<az.n;iz++) {
	    m0[0][ix][iz] = 1.;
	    n0[0][ix][iz] = - c1 * (bb[ix][iz]*bb[ix][iz]) /  aa[ix][iz];
	    r0[0][ix][iz] =   c2 * (bb[ix][iz]*bb[ix][iz]) / (aa[ix][iz] * aa[ix][iz]);

	    m0[0][ix][iz] *= kmu;
	    n0[0][ix][iz] *= knu;
	    r0[0][ix][iz] *= kro;
	}
    }
}

void rweone_ffd_coef(
    float **aa,
    float **bb,
    float **a0,
    float **b0)
/*< FFD coefficients >*/
{
    float ***d1, ***d2;
    int ii,ir,iz,ix;
    float tt,tt2, *tb;
    float boa, boa2, boa4;
    float bob, bob2, bob4;
    float aoa, aoa3;

    d1=sf_floatalloc3(az.n,ax.n,ar.n);
    d2=sf_floatalloc3(az.n,ax.n,ar.n);

    /* find max(b0) */
    tb=sf_floatalloc (az.n*ar.n);    
    ii=0;
    for(ir=0;ir<ar.n;ir++) {
	for(iz=0;iz<az.n;iz++) {
	    tb[ii] = b0[ir][iz];
	    ii++;
	}
    }
    tt = SF_ABS(sf_quantile(ar.n*az.n-1,ar.n*az.n,tb));
    tt2= tt*tt;
    free(tb);

    for(ir=0;ir<ar.n;ir++) {
	for(iz=0;iz<az.n;iz++) {
	    boa = ( b0[ir][iz]/tt ) / a0[ir][iz];
	    boa2= boa *boa ;
	    boa4= boa2*boa2;

	    for(ix=0;ix<ax.n;ix++) {
		bob = bb[ix][iz]/b0[ir][iz];
		bob2= bob*bob;

		aoa = a0[ir][iz]/aa[ix][iz];
		aoa3= aoa*aoa*aoa;

		d1[ir][ix][iz] = a0[ir][iz] * boa2 * ( bob2 * aoa  - 1.);
		d2[ir][ix][iz] = a0[ir][iz] * boa4 * ( bob4 * aoa3 - 1.);
	    }
	}
    }

    /* regularize d1 */
    for(ir=0;ir<ar.n;ir++) {
	for(ix=0;ix<ax.n;ix++) {
	    for(iz=0;iz<az.n;iz++) {
		if( SF_ABS(d1[ir][ix][iz]) < 1e-6) d1[ir][ix][iz]=1e-6;
	    }
	}
    }

    for(ir=0;ir<ar.n;ir++) {
	for(ix=0;ix<ax.n;ix++) {
	    for(iz=0;iz<az.n;iz++) {

		m0[ir][ix][iz] =       d1[ir][ix][iz] / tt2;
		n0[ir][ix][iz] = -c1 * d1[ir][ix][iz] * d1[ir][ix][iz];
		r0[ir][ix][iz] =  c2 * d2[ir][ix][iz];

		m0[ir][ix][iz] *= kmu;
		n0[ir][ix][iz] *= knu;
		r0[ir][ix][iz] *= kro;
	    }
	}
    }
    
    free(**d1); free(*d1); free(d1);
    free(**d2); free(*d2); free(d2);
}

void rweone_psc_coef(
    float **aa,
    float **bb,
    float **a0,
    float **b0)
/*< PSC coefficients >*/
{
    float tt,tt2, *tb;
    int ii,ir,iz,ix;
    float boa, boa2;
    float aoa;
    float bob;

    /* find max(b0) */
    tb=sf_floatalloc (az.n*ar.n);
    ii=0;
    for(ir=0;ir<ar.n;ir++) {
	for(iz=0;iz<az.n;iz++) {
	    tb[ii] = b0[ir][iz];
	    ii++;
	}
    }
    tt = SF_ABS(sf_quantile(ar.n*az.n-1,ar.n*az.n,tb));
    tt2= tt*tt;
    free(tb);

    for(ir=0;ir<ar.n;ir++) {
	for(iz=0;iz<az.n;iz++) {
	    boa = ( b0[ir][iz]/tt ) / a0[ir][iz];
	    boa2=boa*boa;

	    for(ix=0;ix<ax.n;ix++) {
		aoa = aa[ix][iz]/a0[ir][iz];
		bob = bb[ix][iz]/b0[ir][iz];

		m0[ir][ix][iz] = 1 / tt2;
		n0[ir][ix][iz] = a0[ir][iz] * boa2 * (c1*(aoa-1.) - (bob-1.));
		r0[ir][ix][iz] =    3 * c2 * boa2;

		m0[ir][ix][iz] *= kmu;
		n0[ir][ix][iz] *= knu;
		r0[ir][ix][iz] *= kro;
	    }
	}
    }
}

/*------------------------------------------------------------*/
/* 
 * Extrapolation components
 */
/*------------------------------------------------------------*/

void rweone_phs(
    complex float *v,
    float w,
    int iz,
    int ir,
    float **a0,
    float **b0
    )
/*< Fourier-domain phase shift >*/
{
    int ix,ikx;
    float kx,arg;
    float ta,tb,tt;
    complex float ikz;

    rweone_fft(false,v);

    for(ix=0;ix<ax.n;ix++) {
	ikx = KMAP(ix,ax.n);
	kx = okx + ikx * dkx;

	ta = a0[ir][iz]*w;
	tb = b0[ir][iz]*kx;
	tt = tb/ta;
	arg = 1.0 - tt*tt;

	if(arg<0.) {
	    ikz =     ta * sqrtf(-arg);
	} else {
	    ikz = I * ta * sqrtf(+arg);
	}

	v[ix] *= cexp( ikz * (-az.d));
    }
   
    rweone_fft( true,v);
}

void rweone_ssf(
    complex float *v,
    float w,
    int iz,
    int ir,
    float **aa,
    float **a0)
/*< split-step Fourier correction >*/
{
    int ix;
    complex float ikz;

    for(ix=0; ix<ax.n; ix++) {	
	ikz = I * w * (aa[ix][iz] - a0[ir][iz]);
	v[ix] *= cexpf( ikz * (-az.d) );
    }
}

void rweone_ssh(
    complex float *v,
    float w,
    int   iz,
    float **aa)
/*< space-domain phase shift >*/
{
    int ix;

    for(ix=0;ix<ax.n;ix++) {
	v[ix] *= cexpf( -I*w*aa[ix][iz] * az.d);
    }
}

void rweone_fds(
    complex float *v,
    float w,
    int iz,
    int ir)
/*< finite-differences solver >*/
{
    int ix;

    for(ix=0;ix<ax.n;ix++) {
	mu[ix] =  m0[ir][ix][iz];
	nu[ix] = -n0[ir][ix][iz] / w;
	ro[ix] =  r0[ir][ix][iz] /(w*w);
    }

    for(ix=0;ix<ax.n;ix++) {
	lk[ix] = ro[ix]+sixth*mu[ix] + I*nu[ix];
	rk[ix] = ro[ix]+sixth*mu[ix] - I*nu[ix];
	ck[ix] =              mu[ix];
    }

    /*
          b      |         a           |     c
      ro + I nu  | mu - 2(ro + I nu)   | ro + I nu
    */
    for(ix=0;ix<ax.n  ;ix++) { a[ix] = ck[ix] - 2.*lk[ix]; }
    for(ix=0;ix<ax.n-1;ix++) { b[ix] = lk[ix];   }
    for(ix=0;ix<ax.n-1;ix++) { c[ix] = lk[ix+1]; }

    /* 
      ro - I nu  |  mu - 2(ro - I nu)  | ro - I nu
    */
    for(ix=1;ix<ax.n-1;ix++) {
	u[ix] = 
	    rk[ix] *(v[ix-1]+v[ix+1]) +  
	    (ck[ix]-2.0*rk[ix]) * v[ix];
    }
    ix=     0; u[ix]=(ck[ix]-2.0*rk[ix])*v[ix] + rk[ix]*v[ix+1];
    ix=ax.n-1; u[ix]=(ck[ix]-2.0*rk[ix])*v[ix] + rk[ix]*v[ix-1];

    rweone_thr(a,b,c,u,ax.n);

    for(ix=0;ix<ax.n;ix++) { v[ix] = u[ix]; };
}

/*------------------------------------------------------------*/
/* 
 * Utilities
 */
/*------------------------------------------------------------*/

void rweone_thr(
    complex float *a,
    complex float *b,
    complex float *c,
    complex float *v,
    int n)
/*< tridiagonal solver >*/
{
    int i;

    b[0]/=a[0];
    v[0]/=a[0];
    for(i=1;i<n;i++) {
	a[i] -= c[i-1]*b[i-1];
	if(i<n-1) b[i]/=a[i];
	v[i]=( v[i]-c[i-1]*v[i-1] ) / a[i];
    }

    for(i=n-2;i>=0;i--) {
	v[i] -= b[i]*v[i+1];
    }
}

void rweone_tap_init()
/*< initialize taper >*/
{
    int itap,jtap;
    int ix;

    for(ix=0;ix<ax.n;ix++) {
	tap[ix] = 1.;
    }

    if(ntap>1) {
	for(itap=0;itap<ntap;itap++) {
	    jtap = SF_ABS(ntap-1-itap);

	    tap[itap] = cosf(SF_PI/2. * jtap/(ntap-1));
/*	    tap[itap] = sinf(SF_PI/2. * itap/ntap);*/
/*	    tap[itap] *= tap[itap];*/
	    tap[ax.n-1 - itap] = tap[itap];
	}
    }
}

void rweone_tap(complex float *v)
/*< apply taper >*/
{
    int ix;
    for(ix=0;ix<ax.n;ix++) {
	v[ix] *= tap[ix];
    }
}


void rweone_fft_init(int n)
/*< init FFT >*/
{
    if(n>0) {
	ffts = 1./n; 
    } else {
	sf_error("KISS FFT scale error: n==0",__FILE__);
    }

    forw = kiss_fft_alloc(n,0,NULL,NULL);
    invs = kiss_fft_alloc(n,1,NULL,NULL);

    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error",__FILE__);
}

void rweone_fft(
    bool inv,
    complex float *d)
/*< apply FFT >*/
{
    int ix;

    if(inv) for(ix=0;ix<ax.n;ix++) { d[ix] *= ffts; }
    if(inv) {
	kiss_fft(invs,
		 (const kiss_fft_cpx *) d, 
		 (      kiss_fft_cpx *) d);
    } else {
	kiss_fft(forw,
		 (const kiss_fft_cpx *) d, 
		 (      kiss_fft_cpx *) d);
    }
}


void rweone_mrs_init()
/*< init MRS >*/
{
    int ix;
    for(ix=0;ix<ax.n;ix++) {
	msk[ix] = 1.;
	mtt[ix] = 1.;
    }
}

void rweone_mrs(
    complex float *v,
    int iz,
    int ir,
    float **m,
    complex float *d)
/*< combine MRS >*/
{
    int ix, iloop;

    for(ix=0;ix<ax.n;ix++) {
	if(m[ix][iz]-1 == ir) {
	    mtt[ix] = 1.;
	} else {
	    mtt[ix] = 0.;
	}
    }

    for(iloop=0;iloop<nloop;iloop++) {

	msk[0] = mtt[0];
	for(ix=1;ix<ax.n-1;ix++) {
	    msk[ix] = 0.5 * (mtt[ix-1] + mtt[ix+1] );
	}
	msk[ax.n-1] = mtt[ax.n-2];

	for(ix=0;ix<ax.n;ix++) {
	    mtt[ix] = msk[ix];
	}
    }

    for(ix=0;ix<ax.n;ix++) {
	d[ix] += msk[ix] * v[ix];
    }
}

void cwrite(complex float x)
/*< output a complex number >*/
{
    sf_warning("(%f,%f)",crealf(x), cimagf(x));
}
