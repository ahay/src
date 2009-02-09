/* one-way Rienmannian Wavefield Extrapolation 
   extrapolation subroutines
*/
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
/*^*/

#include "rweone.h"

static sf_axa ag,at,aw,ar;
static int method;
static bool verb;

static float ***m0, ***n0, ***r0;
static float   *mu,   *nu,   *ro;

static sf_complex *lk,*rk,*ck;
static sf_complex *a,*b,*c,*u;

static float c1,c2,sixth;
static float kmu,knu,kro;

static float *tap;
static int   ntap;

static kiss_fft_cfg forw, invs;
static float ffts;
static float okg,dkg;

sf_complex *wfl,*wtt;
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
    sf_complex *ddd,
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
	ddd[ig] = sf_cmplx(0.,0.);
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
    sf_complex *ddd,
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
    sf_complex *v,
    float w,
    float a0,
    float b0
    )
/*< Fourier-domain phase shift >*/
{
    int ig,ikg;
    float kg,arg;
    float ta,tb,tt;
    sf_complex ikz;

    rweone_fft(false,(kiss_fft_cpx*) v);

    for(ig=0;ig<ag.n;ig++) {
	ikg = KMAP(ig,ag.n);
	kg = okg + ikg * dkg;
	
	ta = a0*w;
	tb = b0*kg;
	tt = tb/ta;
	arg = 1.0 - tt*tt;
	
	if(arg<0.) {
	    ikz = sf_cmplx(SF_ABS(ta) * sqrtf(-arg),0.);
	} else {
	    ikz = sf_cmplx(0.,ta  * sqrtf(+arg));
	}
#ifdef SF_HAS_COMPLEX_H
	v[ig] *= cexpf( ikz * (-at.d));
#else
	v[ig] = sf_cmul(v[ig],cexpf(sf_crmul(ikz,-at.d)));
#endif
    }
    
    rweone_fft( true,(kiss_fft_cpx*) v);
}

void rweone_phs(
    sf_complex *v,
    float w,
    float a0,
    float b0
    )
/*< Fourier-domain phase shift >*/
{
    int ig,ikg;
    float kg;

    float a2,b2,k2;
    sf_complex iw,ikt,w2;

    a2 = a0*a0;
    b2 = b0*b0;

    iw = sf_cmplx(2e-3,-w);
#ifdef SF_HAS_COMPLEX_H
    w2 = iw*iw;
#else
    w2 = sf_cmul(iw,iw);
#endif

    rweone_fft(false,(kiss_fft_cpx*) v);

    for(ig=0;ig<ag.n;ig++) {
	ikg = KMAP(ig,ag.n);
	kg = okg + ikg * dkg;
	k2 = kg*kg;
	
#ifdef SF_HAS_COMPLEX_H
	ikt = csqrtf( w2*a2 + k2*b2 );
	v[ig] *= cexpf(-ikt*at.d);
#else
	ikt = csqrtf(sf_cadd(sf_crmul(w2,a2),sf_cmplx(k2*b2,0.)));
	v[ig] = sf_cmul(v[ig],cexpf(sf_crmul(ikt,-at.d)));
#endif
    }
    
    rweone_fft( true,(kiss_fft_cpx*) v);
}


void rweone_ssf(
    sf_complex *v,
    float w,
    float *aa,
    float  a0)
/*< split-step Fourier correction >*/
{
    int ig;
    sf_complex ikz;

    for(ig=0; ig<ag.n; ig++) {	
	ikz = sf_cmplx(0.,w * (aa[ig] - a0));
#ifdef SF_HAS_COMPLEX_H
	v[ig] *= cexpf( ikz * ( at.d) );
#else
	v[ig] = sf_cmul(v[ig],cexpf(sf_crmul(ikz, at.d)));
#endif
    }
}

void rweone_ssh(
    sf_complex *v,
    float w,
    float *aa)
/*< space-domain phase shift >*/
{
    int ig;
    sf_complex ikz;

    for(ig=0;ig<ag.n;ig++) {
	ikz = sf_cmplx(0.,w * aa[ig]);
#ifdef SF_HAS_COMPLEX_H
	v[ig] *= cexpf( ikz * at.d );
#else
	v[ig] = sf_cmul(v[ig],cexpf(sf_crmul( ikz, at.d )));
#endif
    }
}

void rweone_fds(
    sf_complex *v,
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
	lk[ig] = sf_cmplx(ro[ig]+sixth*mu[ig],+nu[ig]);
	rk[ig] = sf_cmplx(ro[ig]+sixth*mu[ig],-nu[ig]);
	ck[ig] = sf_cmplx(             mu[ig],0.);
    }

    /*
          b      |         a           |     c
      ro + I nu  | mu - 2(ro + I nu)   | ro + I nu
    */
#ifdef SF_HAS_COMPLEX_H
    for(ig=0;ig<ag.n  ;ig++) { a[ig] = ck[ig] - 2.*lk[ig]; }
#else
    for(ig=0;ig<ag.n  ;ig++) { a[ig] = sf_csub(ck[ig],sf_crmul(lk[ig],2.)); } 
#endif
    for(ig=0;ig<ag.n-1;ig++) { b[ig] = lk[ig];   }
    for(ig=0;ig<ag.n-1;ig++) { c[ig] = lk[ig+1]; }

    /* 
      ro - I nu  |  mu - 2(ro - I nu)  | ro - I nu
    */
    for(ig=1;ig<ag.n-1;ig++) {
#ifdef SF_HAS_COMPLEX_H
	u[ig] = 
	    rk[ig] *(v[ig-1]+v[ig+1]) +  
	    (ck[ig]-2.0*rk[ig]) * v[ig];
#else
	u[ig] = 
	    sf_cadd(sf_cmul(rk[ig],sf_cadd(v[ig-1],v[ig+1])),  
		    sf_cmul(sf_csub(ck[ig],sf_crmul(rk[ig],2.)),v[ig]));
#endif	
    }
#ifdef SF_HAS_COMPLEX_H
    ig=     0; u[ig]=(ck[ig]-2.0*rk[ig])*v[ig] + rk[ig]*v[ig+1];
    ig=ag.n-1; u[ig]=(ck[ig]-2.0*rk[ig])*v[ig] + rk[ig]*v[ig-1];
#else
    ig=     0; u[ig]=sf_cadd(
	sf_cmul(sf_csub(ck[ig],sf_crmul(rk[ig],2.0)),v[ig]),
	sf_cmul(rk[ig],v[ig+1]));
    ig=ag.n-1; u[ig]=sf_cadd(
	sf_cmul(sf_csub(ck[ig],sf_crmul(rk[ig],2.0)),v[ig]),
	sf_cmul(rk[ig],v[ig-1]));
#endif
    
    rweone_thr(a,b,c,u,ag.n);

    for(ig=0;ig<ag.n;ig++) { v[ig] = u[ig]; };
}

/*------------------------------------------------------------*/
/* 
 * Utilities
 */
/*------------------------------------------------------------*/

void rweone_thr(
    sf_complex *a,
    sf_complex *b,
    sf_complex *c,
    sf_complex *v,
    int n)
/*< tridiagonal solver >*/
{
    int i;

#ifdef SF_HAS_COMPLEX_H
    b[0]/=a[0];
    v[0]/=a[0];
#else
    b[0] = sf_cdiv(b[0],a[0]);
    v[0] = sf_cdiv(v[0],a[0]);
#endif
    for(i=1;i<n;i++) {
#ifdef SF_HAS_COMPLEX_H
	a[i] -= c[i-1]*b[i-1];
	if(i<n-1) b[i]/=a[i];
	v[i]=( v[i]-c[i-1]*v[i-1] ) / a[i];
#else
	a[i] = sf_csub(a[i],sf_cmul(c[i-1],b[i-1]));
	if(i<n-1) b[i] = sf_cdiv(b[i],a[i]);
	v[i]=sf_cdiv(sf_csub(v[i],sf_cmul(c[i-1],v[i-1])),a[i]);
#endif
    }

    for(i=n-2;i>=0;i--) {
#ifdef SF_HAS_COMPLEX_H
	v[i] -= b[i]*v[i+1];
#else
	v[i] = sf_csub(v[i],sf_cmul(b[i],v[i+1]));
#endif
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

void rweone_tap(sf_complex *v)
/*< apply taper >*/
{
    int ig;
    for(ig=0;ig<ag.n;ig++) {
#ifdef SF_HAS_COMPLEX_H
	v[ig] *= tap[ig];
#else
	v[ig] = sf_crmul(v[ig],tap[ig]);
#endif
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
    kiss_fft_cpx *d)
/*< apply FFT >*/
{
    int ig;

    if(inv) {

	kiss_fft(invs,d,d);
	for(ig=0;ig<ag.n;ig++) { d[ig] = sf_crmul(d[ig],ffts); }

    } else {

	for(ig=0;ig<ag.n;ig++) { d[ig] = sf_crmul(d[ig],ffts); }
	kiss_fft(forw,d,d);

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
    sf_complex *v,
    sf_complex *d,
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
#ifdef SF_HAS_COMPLEX_H
	d[ig] += msk[ig] * v[ig];
#else
	d[ig] = sf_cadd(d[ig],sf_crmul(v[ig],msk[ig]));
#endif
    }
}

/*------------------------------------------------------------*/

void rweone_zoi(
    bool        adj,
    sf_complex *ddd,
    float      *iii)
/*< zero-offset imaging condition >*/
{
    int ig;

    if(adj) { /* modeling */
	for(ig=0;ig<ag.n;ig++) {
#ifdef SF_HAS_COMPLEX_H
	    ddd[ig]   += iii[ig];
#else
	    ddd[ig].r += iii[ig];
#endif
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
    sf_complex *swf,
    sf_complex *rwf,
    float         *iii)
/*< shot-record imaging condition >*/
{
    int ig;

    rweone_tap(swf);
    rweone_tap(rwf);
    for(ig=0;ig<ag.n;ig++) {
#ifdef SF_HAS_COMPLEX_H
	iii[ig] += crealf(        conjf(swf[ig]) * rwf[ig]  );
#else
	iii[ig] += crealf(sf_cmul(conjf(swf[ig]),  rwf[ig]) );
#endif
    }
}
