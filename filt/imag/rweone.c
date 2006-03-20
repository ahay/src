/* one-way Rienmannian Wavefield Extrapolation 
   extrapolation subroutines
*/
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

#include "rweone.h"

static sf_axa ag,at,aw,ar;
static int method;
static bool verb;

static float ***m0, ***n0, ***r0;
static float   *mu,   *nu,   *ro;

static complex float *lk,*rk,*ck;
static complex float *a,*b,*c,*u;

static float c1,c2,sixth;
static float kmu,knu,kro;

static float *tap;
static int   ntap;

static kiss_fft_cfg forw, invs;
static float ffts;
static float okg,dkg;

complex float *wfl,*wtt;
float *mtt,*msk;
int nloop;

#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);

void rweone_init(
    sf_axis ag_,
    sf_axis at_,
    sf_axis aw_,
    sf_axis ar_,
    int method_,
    bool verb_)
/*< initialize >*/
{
    ag = sf_nod(ag_);
    at = sf_nod(at_);
    aw = sf_nod(aw_);
    ar = sf_nod(ar_);

    method = method_;
    verb   = verb_;

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    /* space for F-D */
    m0 = sf_floatalloc3(ag.n,ar.n,at.n);
    n0 = sf_floatalloc3(ag.n,ar.n,at.n);
    r0 = sf_floatalloc3(ag.n,ar.n,at.n);
    
    mu = sf_floatalloc(ag.n);
    nu = sf_floatalloc(ag.n);
    ro = sf_floatalloc(ag.n);

    if(! sf_getfloat("c1",&c1)) c1=0.47;
    if(! sf_getfloat("c2",&c2)) c2=0.37;
    if(! sf_getfloat("sixth",&sixth)) sixth=0.122996;

    lk = sf_complexalloc(ag.n);
    rk = sf_complexalloc(ag.n);
    ck = sf_complexalloc(ag.n);

    a = sf_complexalloc(ag.n);
    b = sf_complexalloc(ag.n-1);
    c = sf_complexalloc(ag.n-1);
    u = sf_complexalloc(ag.n);

    kmu = 1./ at.d;             /* 1/ dz */
    knu = 1 /(   2*ag.d*ag.d);  /* 1/( 2 dx^2) */
    kro = 1./(at.d*ag.d*ag.d);  /* 1/(dz dx^2) */

    /* taper */
    if(! sf_getint("ntap",&ntap)) ntap=1;
    tap = sf_floatalloc(ag.n);
    rweone_tap_init();

    /* Fourier-domain */
    rweone_fft_init(ag.n);
    okg = -  SF_PI /        ag.d ;
    dkg = 2.*SF_PI /(ag.n * ag.d);

    /* MRS temp arrays */
    wfl = sf_complexalloc(ag.n);
    wtt = sf_complexalloc(ag.n);

    msk = sf_floatalloc(ag.n);
    mtt = sf_floatalloc(ag.n);

    if(! sf_getint("nloop",&nloop)) nloop=3;
    rweone_mrs_init();
}

void rweone_fk(
    float w,
    complex float *ddd,
    float *aa,
    float *a0,
    float *b0,
    float *mm,
    int it)
/*< one F-K extrapolation step >*/
{
    int ig,ir;

    for(ig=0;ig<ag.n;ig++) {
	wtt[ig] = ddd[ig];
	ddd[ig] = 0.;
    }
    for(ir=0;ir<ar.n;ir++) {
	for(ig=0;ig<ag.n;ig++) {
	    wfl[ig] = wtt[ig];
	}
	
	rweone_phs(wfl,w,a0[ir],b0[ir]);
	rweone_ssf(wfl,w,aa,a0[ir]);
	/* skip F-D for SSF */
	if(method!=1) { rweone_fds(wfl,w,it,ir); } 
	rweone_mrs(wfl,ddd,mm,ir);
    }
    
}

void rweone_fx(
    float w,
    complex float *ddd,
    float *aa,
    int it)
/*< one F-X extrapolation step >*/
{
    rweone_ssh(ddd,w,aa);     /* space shift */
    rweone_fds(ddd,w,it, 0);  /* finite-differences */
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
    int it,ig;

    if(verb) sf_warning("compute XFD coef");

    for(it=0;it<at.n;it++) {
	for(ig=0;ig<ag.n;ig++) {
	    m0[it][0][ig] = 1.;
	    n0[it][0][ig] = - c1 * (bb[it][ig]*bb[it][ig]) /  aa[it][ig];
	    r0[it][0][ig] =   c2 * (bb[it][ig]*bb[it][ig]) / (aa[it][ig] * aa[it][ig]);

	    m0[it][0][ig] *= kmu;
	    n0[it][0][ig] *= knu;
	    r0[it][0][ig] *= kro;
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
    int ii,ir,it,ig;
    float tt,tt2, *tb;
    float boa, boa2, boa4;
    float bob, bob2, bob4;
    float aoa, aoa3;

    if(verb) sf_warning("compute FFD coef");

    d1=sf_floatalloc3(ag.n,ar.n,at.n);
    d2=sf_floatalloc3(ag.n,ar.n,at.n);

    /* find max(b0) */
    tb=sf_floatalloc (ar.n*at.n); 
    ii=0;
    for(it=0;it<at.n;it++) {
	for(ir=0;ir<ar.n;ir++) {
	    tb[ii] = b0[it][ir];
	    ii++;
	}
    }
    tt = SF_ABS(sf_quantile(ar.n*at.n-1,ar.n*at.n,tb));
    tt2= tt*tt;
    free(tb);

    for(it=0;it<at.n;it++) {
	for(ir=0;ir<ar.n;ir++) {
	    boa = ( b0[it][ir]/tt ) / a0[it][ir];
	    boa2= boa *boa ;
	    boa4= boa2*boa2;

	    for(ig=0;ig<ag.n;ig++) {
		bob = bb[it][ig]/b0[it][ir];
		bob2= bob *bob;
		bob4= bob2*bob2;

		aoa = a0[it][ir]/aa[it][ig];
		aoa3= aoa*aoa*aoa;

		d1[it][ir][ig] = a0[it][ir] * boa2 * ( bob2 * aoa  - 1.);
		d2[it][ir][ig] = a0[it][ir] * boa4 * ( bob4 * aoa3 - 1.);
	    }
	}
    }

    /* regularize d1 */
    for(it=0;it<at.n;it++) {
	for(ir=0;ir<ar.n;ir++) {
	    for(ig=0;ig<ag.n;ig++) {
		if( SF_ABS(d1[it][ir][ig]) < 1e-6) d1[it][ir][ig]=1e-6;
	    }
	}
    }
    
    for(it=0;it<at.n;it++) {
	for(ir=0;ir<ar.n;ir++) {
	    for(ig=0;ig<ag.n;ig++) {

		m0[it][ir][ig] =       d1[it][ir][ig] / tt2;
		n0[it][ir][ig] = -c1 * d1[it][ir][ig] * d1[it][ir][ig];
		r0[it][ir][ig] =  c2 * d2[it][ir][ig];

		m0[it][ir][ig] *= kmu;
		n0[it][ir][ig] *= knu;
		r0[it][ir][ig] *= kro;
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
    int ii,ir,it,ig;
    float boa, boa2;
    float aoa;
    float bob;

    if(verb) sf_warning("compute PSC coef");

    /* find max(b0) */
    tb=sf_floatalloc (ar.n*at.n);
    ii=0;
    for(it=0;it<at.n;it++) {
	for(ir=0;ir<ar.n;ir++) {
	    tb[ii] = b0[it][ir];
	    ii++;
	}
    }
    tt = SF_ABS(sf_quantile(ar.n*at.n-1,ar.n*at.n,tb));
    tt2= tt*tt;
    free(tb);

    for(it=0;it<at.n;it++) {
	for(ir=0;ir<ar.n;ir++) {
	    boa = ( b0[it][ir]/tt ) / a0[it][ir];
	    boa2=boa*boa;

	    for(ig=0;ig<ag.n;ig++) {
		aoa = aa[it][ig]/a0[it][ir];
		bob = bb[it][ig]/b0[it][ir];

		m0[it][ir][ig] = 1 / tt2;
		n0[it][ir][ig] = a0[it][ir] * boa2 * (c1*(aoa-1.) - (bob-1.));
		r0[it][ir][ig] =    3 * c2 * boa2;

		m0[it][ir][ig] *= kmu;
		n0[it][ir][ig] *= knu;
		r0[it][ir][ig] *= kro;
	    }
	}
    }
}

/*------------------------------------------------------------*/
/* 
 * Extrapolation components
 */
/*------------------------------------------------------------*/

void rweone_phs_old(
    complex float *v,
    float w,
    float a0,
    float b0
    )
/*< Fourier-domain phase shift >*/
{
    int ig,ikg;
    float kg,arg;
    float ta,tb,tt;
    complex float ikz;

    rweone_fft(false,v);

    for(ig=0;ig<ag.n;ig++) {
	ikg = KMAP(ig,ag.n);
	kg = okg + ikg * dkg;
	
	ta = a0*w;
	tb = b0*kg;
	tt = tb/ta;
	arg = 1.0 - tt*tt;
	
	if(arg<0.) {
	    ikz = SF_ABS(ta) * sqrtf(-arg);
	} else {
	    ikz =    I * ta  * sqrtf(+arg);
	}
	v[ig] *= cexp( ikz * (-at.d));
    }
    
    rweone_fft( true,v);
}

void rweone_phs(
    complex float *v,
    float w,
    float a0,
    float b0
    )
/*< Fourier-domain phase shift >*/
{
    int ig,ikg;
    float kg;

    float a2,b2,k2;
    complex float iw,ikt,w2;

    a2 = a0*a0;
    b2 = b0*b0;

    iw = 2e-3 - I * w;
    w2 = iw*iw;

    rweone_fft(false,v);

    for(ig=0;ig<ag.n;ig++) {
	ikg = KMAP(ig,ag.n);
	kg = okg + ikg * dkg;
	k2 = kg*kg;
	
	ikt = csqrtf( w2*a2 + k2*b2 );

	v[ig] *= cexp(-ikt*at.d);
    }
    
    rweone_fft( true,v);
}


void rweone_ssf(
    complex float *v,
    float w,
    float *aa,
    float  a0)
/*< split-step Fourier correction >*/
{
    int ig;
    complex float ikz;

    for(ig=0; ig<ag.n; ig++) {	
	ikz = I * w * (aa[ig] - a0);
	v[ig] *= cexpf( ikz * (-at.d) );
    }
}

void rweone_ssh(
    complex float *v,
    float w,
    float *aa)
/*< space-domain phase shift >*/
{
    int ig;
    complex float ikz;

    for(ig=0;ig<ag.n;ig++) {
	ikz = I * w * aa[ig];
	v[ig] *= cexpf( ikz * at.d );
    }
}

void rweone_fds(
    complex float *v,
    float w,
    int it,
    int ir)
/*< finite-differences solver >*/
{
    int ig;

    w *=-1; /* fix sign convention */

    for(ig=0;ig<ag.n;ig++) {
	mu[ig] =  m0[it][ir][ig];
	nu[ig] = -n0[it][ir][ig] / w;
	ro[ig] =  r0[it][ir][ig] /(w*w);
    }

    for(ig=0;ig<ag.n;ig++) {
	lk[ig] = ro[ig]+sixth*mu[ig] + I*nu[ig];
	rk[ig] = ro[ig]+sixth*mu[ig] - I*nu[ig];
	ck[ig] =              mu[ig];
    }

    /*
          b      |         a           |     c
      ro + I nu  | mu - 2(ro + I nu)   | ro + I nu
    */
    for(ig=0;ig<ag.n  ;ig++) { a[ig] = ck[ig] - 2.*lk[ig]; }
    for(ig=0;ig<ag.n-1;ig++) { b[ig] = lk[ig];   }
    for(ig=0;ig<ag.n-1;ig++) { c[ig] = lk[ig+1]; }

    /* 
      ro - I nu  |  mu - 2(ro - I nu)  | ro - I nu
    */
    for(ig=1;ig<ag.n-1;ig++) {
	u[ig] = 
	    rk[ig] *(v[ig-1]+v[ig+1]) +  
	    (ck[ig]-2.0*rk[ig]) * v[ig];
    }
    ig=     0; u[ig]=(ck[ig]-2.0*rk[ig])*v[ig] + rk[ig]*v[ig+1];
    ig=ag.n-1; u[ig]=(ck[ig]-2.0*rk[ig])*v[ig] + rk[ig]*v[ig-1];

    rweone_thr(a,b,c,u,ag.n);

    for(ig=0;ig<ag.n;ig++) { v[ig] = u[ig]; };
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
    int ig;

    for(ig=0;ig<ag.n;ig++) {
	tap[ig] = 1.;
    }

    if(ntap>1) {
	for(itap=0;itap<ntap;itap++) {
	    jtap = SF_ABS(ntap-1-itap);

	    tap[itap] = cosf(SF_PI/2. * jtap/(ntap-1));
/*	    tap[itap] = sinf(SF_PI/2. * itap/ntap);*/
/*	    tap[itap] *= tap[itap];*/
	    tap[ag.n-1 - itap] = tap[itap];
	}
    }
}

void rweone_tap(complex float *v)
/*< apply taper >*/
{
    int ig;
    for(ig=0;ig<ag.n;ig++) {
	v[ig] *= tap[ig];
    }
}


void rweone_fft_init(int n)
/*< init FFT >*/
{
    if(n>0) {
	ffts = 1./sqrt(n); 
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
    int ig;

    if(inv) {

	kiss_fft(invs,
		 (const kiss_fft_cpx *) d, 
		 (      kiss_fft_cpx *) d);
	for(ig=0;ig<ag.n;ig++) { d[ig] *= ffts; }

    } else {

	for(ig=0;ig<ag.n;ig++) { d[ig] *= ffts; }
	kiss_fft(forw,
		 (const kiss_fft_cpx *) d, 
		 (      kiss_fft_cpx *) d);
    }
}


void rweone_mrs_init()
/*< init MRS >*/
{
    int ig;
    for(ig=0;ig<ag.n;ig++) {
	msk[ig] = 1.;
	mtt[ig] = 1.;
    }
}

void rweone_mrs(
    complex float *v,
    complex float *d,
    float *m,
    int ir)
/*< combine MRS >*/
{
    int ig, iloop;

    for(ig=0;ig<ag.n;ig++) {
	if(m[ig]-1 == ir) {
	    mtt[ig] = 1.;
	} else {
	    mtt[ig] = 0.;
	}
    }

    for(iloop=0;iloop<nloop;iloop++) {
	msk[0] = mtt[1];
	for(ig=1;ig<ag.n-1;ig++) {
	    msk[ig] = 0.5 * (mtt[ig-1] + mtt[ig+1] );
	}
	msk[ag.n-1] = mtt[ag.n-2];

	for(ig=0;ig<ag.n;ig++) {
	    mtt[ig] = msk[ig];
	}
    }

    for(ig=0;ig<ag.n;ig++) {
	d[ig] += msk[ig] * v[ig];
    }
}

/*------------------------------------------------------------*/

void cwrite(complex float x)
/*< output a complex number >*/
{
    sf_warning("(%f,%f)",crealf(x), cimagf(x));
}

/*------------------------------------------------------------*/

void rweone_zoi(
    bool           adj,
    complex float *ddd,
    float         *iii)
/*< zero-offset imaging condition >*/
{
    int ig;

    if(adj) { /* modeling */
	for(ig=0;ig<ag.n;ig++) {
	    ddd[ig] += iii[ig];
	}
	rweone_tap(ddd);
    } else {
	rweone_tap(ddd);
	for(ig=0;ig<ag.n;ig++) {
	    iii[ig] += crealf(ddd[ig]);
	}
    }
}

void rweone_spi(
    complex float *swf,
    complex float *rwf,
    float         *iii)
/*< shot-record imaging condition >*/
{
    int ig;

    rweone_tap(swf);
    rweone_tap(rwf);
    for(ig=0;ig<ag.n;ig++) {
	iii[ig] += crealf( conjf(swf[ig]) * rwf[ig] );
    }
}
