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

#include "rweone.h"

static axa ag,at,aw,ar;
static int method;

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
    axa ag_,
    axa at_,
    axa aw_,
    axa ar_,
    int method_)
/*< initialize >*/
{
    ag = ag_;
    at = at_;
    aw = aw_;
    ar = ar_;

    method = method_;
    
    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;
    
    /* space for F-D */
    m0 = sf_floatalloc3(at.n,ag.n,ar.n);
    n0 = sf_floatalloc3(at.n,ag.n,ar.n);
    r0 = sf_floatalloc3(at.n,ag.n,ar.n);
    
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

void rweone_main(
    bool            inv,
    complex float **dat,
    float         **img,
    float         ** aa,
    float         ** bb,
    float         ** mm,
    float         ** a0,
    float         ** b0)
/*< run one-way extrapolation >*/
{
    int iw,it;
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
		if(method==3)      { sf_warning("PSC %d %d",iw,aw.n); }
		else if(method==2) { sf_warning("FFD %d %d",iw,aw.n); }
		else               { sf_warning("SSF %d %d",iw,aw.n); }

		if(inv) { /* modeling */
		    w=aw.o+iw*aw.d; w*=2;
		    w*=-1; /* causal */

		    for(it=at.n-1;it>=0;it--) {
			rweone_fk(w,dat[iw],it,aa,a0,b0,mm);
			rweone_img(inv,dat[iw],img[it]);
		    }
		} else {  /* migration */
		    w=aw.o+iw*aw.d; w*=2;
		    w*=+1; /* anti-causal */

		    for(it=0;it<at.n;it++) {
			rweone_img(inv,dat[iw],img[it]);
			rweone_fk(w,dat[iw],it,aa,a0,b0,mm);
		    }
		}
	    }
	    break;
	case 0:
	default:
	    for(iw=0;iw<aw.n;iw++) {
		sf_warning("XFD %d %d",iw,aw.n);
		
		if(inv) { /* modeling */
		    w=aw.o+iw*aw.d; w*=2;    
		    w*=-1; /* causal */

		    for(it=at.n-1;it>=0;it--) {
			rweone_fx(w,dat[iw],it,aa);
			rweone_img(inv,dat[iw],img[it]);
		    }
		} else { /* migration */
		    w=aw.o+iw*aw.d; w*=2;    
		    w*=+1; /* anti-causal */

		    for(it=0;it<at.n;it++) {
			rweone_img(inv,dat[iw],img[it]);
			rweone_fx(w,dat[iw],it,aa);
		    }
		}
	    }
	    break;
    } /* end switch method */
}

/*------------------------------------------------------------*/

void rweone_img(
    bool           inv,
    complex float *ddd,
    float         *iii)
/*< imaging condition >*/
{
    int ig;

    if(inv) { /* modeling */
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

void rweone_fk(
    float w,
    complex float *ddd,
    int it,
    float **aa,
    float **a0,
    float **b0,
    float **mm)
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
	
	rweone_phs(wfl,w,it,ir,a0,b0);
	rweone_ssf(wfl,w,it,ir,aa,a0);
	/* skip F-D for SSF */
	if(method!=1) { rweone_fds(wfl,w,it,ir); } 
	rweone_mrs(wfl,  it,ir,mm,ddd);
    }
    
}

void rweone_fx(
    float w,
    complex float *ddd,
    int it,
    float **aa)
/*< one F-X extrapolation step >*/
{
    rweone_ssh(ddd,w,it,aa);
    rweone_fds(ddd,w,it, 0);
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

    for(ig=0;ig<ag.n;ig++) {
	for(it=0;it<at.n;it++) {
	    m0[0][ig][it] = 1.;
	    n0[0][ig][it] = - c1 * (bb[ig][it]*bb[ig][it]) /  aa[ig][it];
	    r0[0][ig][it] =   c2 * (bb[ig][it]*bb[ig][it]) / (aa[ig][it] * aa[ig][it]);

	    m0[0][ig][it] *= kmu;
	    n0[0][ig][it] *= knu;
	    r0[0][ig][it] *= kro;
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

    d1=sf_floatalloc3(at.n,ag.n,ar.n);
    d2=sf_floatalloc3(at.n,ag.n,ar.n);

    /* find max(b0) */
    tb=sf_floatalloc (at.n*ar.n);    
    ii=0;
    for(ir=0;ir<ar.n;ir++) {
	for(it=0;it<at.n;it++) {
	    tb[ii] = b0[ir][it];
	    ii++;
	}
    }
    tt = SF_ABS(sf_quantile(ar.n*at.n-1,ar.n*at.n,tb));
    tt2= tt*tt;
    free(tb);

    for(ir=0;ir<ar.n;ir++) {
	for(it=0;it<at.n;it++) {
	    boa = ( b0[ir][it]/tt ) / a0[ir][it];
	    boa2= boa *boa ;
	    boa4= boa2*boa2;

	    for(ig=0;ig<ag.n;ig++) {
		bob = bb[ig][it]/b0[ir][it];
		bob2= bob *bob;
		bob4= bob2*bob2;

		aoa = a0[ir][it]/aa[ig][it];
		aoa3= aoa*aoa*aoa;

		d1[ir][ig][it] = a0[ir][it] * boa2 * ( bob2 * aoa  - 1.);
		d2[ir][ig][it] = a0[ir][it] * boa4 * ( bob4 * aoa3 - 1.);
	    }
	}
    }

    /* regularize d1 */
    for(ir=0;ir<ar.n;ir++) {
	for(ig=0;ig<ag.n;ig++) {
	    for(it=0;it<at.n;it++) {
		if( SF_ABS(d1[ir][ig][it]) < 1e-6) d1[ir][ig][it]=1e-6;
	    }
	}
    }

    for(ir=0;ir<ar.n;ir++) {
	for(ig=0;ig<ag.n;ig++) {
	    for(it=0;it<at.n;it++) {

		m0[ir][ig][it] =       d1[ir][ig][it] / tt2;
		n0[ir][ig][it] = -c1 * d1[ir][ig][it] * d1[ir][ig][it];
		r0[ir][ig][it] =  c2 * d2[ir][ig][it];

		m0[ir][ig][it] *= kmu;
		n0[ir][ig][it] *= knu;
		r0[ir][ig][it] *= kro;
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

    /* find max(b0) */
    tb=sf_floatalloc (at.n*ar.n);
    ii=0;
    for(ir=0;ir<ar.n;ir++) {
	for(it=0;it<at.n;it++) {
	    tb[ii] = b0[ir][it];
	    ii++;
	}
    }
    tt = SF_ABS(sf_quantile(ar.n*at.n-1,ar.n*at.n,tb));
    tt2= tt*tt;
    free(tb);

    for(ir=0;ir<ar.n;ir++) {
	for(it=0;it<at.n;it++) {
	    boa = ( b0[ir][it]/tt ) / a0[ir][it];
	    boa2=boa*boa;

	    for(ig=0;ig<ag.n;ig++) {
		aoa = aa[ig][it]/a0[ir][it];
		bob = bb[ig][it]/b0[ir][it];

		m0[ir][ig][it] = 1 / tt2;
		n0[ir][ig][it] = a0[ir][it] * boa2 * (c1*(aoa-1.) - (bob-1.));
		r0[ir][ig][it] =    3 * c2 * boa2;

		m0[ir][ig][it] *= kmu;
		n0[ir][ig][it] *= knu;
		r0[ir][ig][it] *= kro;
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
    int it,
    int ir,
    float **a0,
    float **b0
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

	ta = a0[ir][it]*w;
	tb = b0[ir][it]*kg;
	tt = tb/ta;
	arg = 1.0 - tt*tt;

	if(arg<0.) {
	    ikz =     ta * sqrtf(-arg);
	} else {
	    ikz = I * ta * sqrtf(+arg);
	}

	v[ig] *= cexp( ikz * (-at.d));
    }
   
    rweone_fft( true,v);
}

void rweone_ssf(
    complex float *v,
    float w,
    int it,
    int ir,
    float **aa,
    float **a0)
/*< split-step Fourier correction >*/
{
    int ig;
    complex float ikz;

    for(ig=0; ig<ag.n; ig++) {	
	ikz = I * w * (aa[ig][it] - a0[ir][it]);
	v[ig] *= cexpf( ikz * (-at.d) );
    }
}

void rweone_ssh(
    complex float *v,
    float w,
    int   it,
    float **aa)
/*< space-domain phase shift >*/
{
    int ig;
    complex float ikz;

    for(ig=0;ig<ag.n;ig++) {
	ikz = I * w * aa[ig][it];
	v[ig] *= cexpf( ikz * (-at.d) );
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

    for(ig=0;ig<ag.n;ig++) {
	mu[ig] =  m0[ir][ig][it];
	nu[ig] = -n0[ir][ig][it] / w;
	ro[ig] =  r0[ir][ig][it] /(w*w);
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
    int ig;

    if(inv) for(ig=0;ig<ag.n;ig++) { d[ig] *= ffts; }
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
    int ig;
    for(ig=0;ig<ag.n;ig++) {
	msk[ig] = 1.;
	mtt[ig] = 1.;
    }
}

void rweone_mrs(
    complex float *v,
    int it,
    int ir,
    float **m,
    complex float *d)
/*< combine MRS >*/
{
    int ig, iloop;

    for(ig=0;ig<ag.n;ig++) {
	if(m[ig][it]-1 == ir) {
	    mtt[ig] = 1.;
	} else {
	    mtt[ig] = 0.;
	}
    }

    for(iloop=0;iloop<nloop;iloop++) {

	msk[0] = mtt[0];
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

void cwrite(complex float x)
/*< output a complex number >*/
{
    sf_warning("(%f,%f)",crealf(x), cimagf(x));
}
