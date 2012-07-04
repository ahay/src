/* Helmholtz factorization */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "xkolmog.h"
#include "chelix.h"
#include "cpolydiv.h"
#include "chelicon.h"

static void init_wave(int init, 
		      int nx, float dx,
		      int nz, float dz,
		      sf_complex *pp /* [nx] */,
		      float wov, int nw, int iw)
{
    int ix;
    float x,x0,z0,phase,amp;
    
    x0 = nx*dx/3;
    z0 = nz*dz/3;

    switch(init) {
	case 1: /*  planar wave @ 15deg */
	    for (ix=0; ix < nx; ix++) {
		x = (ix+1)*dx - x0;
		phase  = wov*x*sinf(15*SF_PI/180.);
		pp[ix] = cexpf(sf_cmplx(0.,phase));
	    }
	    break;
	case 2: /* expanding spherical wave */
	    for (ix=0; ix < nx; ix++) {
		x = (ix+1)*dx - x0;
		phase  = wov*hypotf(z0,x);
		pp[ix] = cexpf(sf_cmplx(0.,phase));
	    }
	    break;
	case 3: /* point source */
	    for (ix=0; ix < nx; ix++) {
		pp[ix]=sf_cmplx(0.,0.);
	    }
	    pp[nx/3-1] = sf_cmplx(1.,0.);
	    break;
	case 4: /* collapsing spherical wave */
	    for (ix=0; ix < nx; ix++) {
		x = (ix+1)*dx - x0;
		phase  = -wov*hypotf(z0,x);
		pp[ix] = cexpf(sf_cmplx(0.,phase));
	    }
	    break;
	default:
	    sf_error("Unknown init=%d",init);
    }

    amp = (nw-iw+1.0)/nw;
    amp = cosf((1-amp)*(0.5*SF_PI));
    amp *= amp;

    for (ix=0; ix < nx; ix++) {
#ifdef SF_HAS_COMPLEX_H
	pp[ix] *= amp;
#else
	pp[ix] = sf_crmul(pp[ix],amp);
#endif
    }
}
  
static int pad2(int n)
{
    int p;
    p = 1;
    while ( p < n )
	p *=2;
    return p;
}

int main(int argc, char* argv[])
{
    bool impresp;
    int nt,nx,nz,nw,init,i,padfactor,nfilt,nkol,it,ix,iz,iw;
    float v,dx,dz,lambda,sixth,gamma,epsdamp,pi2,dw,dt, w,wov;
    sf_complex wov2, a, b, c, d, cshift;
    float ***ppp;
    sf_complex *pp, *qq;
    cfilter aa, fac1, fac2;
    sf_file out, imp=NULL;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");

    if (!sf_getint("nz",&nz)) nz=96;
    if (!sf_getint("nx",&nx)) nx=48;
    if (!sf_getint("nt",&nt)) nt=12;
    if (!sf_getint("nw",&nw)) nw=2;
    if (!sf_getint("init",&init)) init=1;

    if (!sf_getfloat("v",&v)) v=1.;
    if (!sf_getfloat("dz",&dz)) dz=1.;
    if (!sf_getfloat("dx",&dx)) dx=2.;
    if (!sf_getfloat("lambda",&lambda)) lambda=nz*dz/4.;

    sf_putint(out,"n1",nz);
    sf_putint(out,"n2",nx);
    sf_putint(out,"n3",nt);

    aa = allocatechelix(9);

    if (!sf_getfloat("sixth",&sixth)) sixth=0.0833;
    if (!sf_getfloat("gamma",&gamma)) gamma=0.667;
    if (!sf_getfloat("epsdamp",&epsdamp)) epsdamp=0.01;
    if (!sf_getint("padfactor",&padfactor)) padfactor=1024;
    if (!sf_getint("nfilt",&nfilt)) nfilt=nx+2;
    if (!sf_getbool("impresp",&impresp)) impresp=false;

    if (impresp) {
	imp = sf_output("imp");
	sf_setformat(imp,"native_complex");
	sf_putint(imp,"n1",2*nx);
	sf_putint(imp,"n2",2);
    }

    ppp = sf_floatalloc3(nz,nx,nt);
    pp = sf_complexalloc(nx*nz);
    qq = sf_complexalloc(nx*nz);

    pi2 = 2.*SF_PI;
    dw = v*pi2/lambda;
    dt =   pi2/(nt*dw);

    nkol=pad2(padfactor*nx);
    /* dkol=pi2/nkol; */

    for (it=0; it < nt; it++) {
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		ppp[it][ix][iz] = 0.;
	    }
	}
    }

    fac1 = allocatechelix(nfilt);
    fac2 = allocatechelix(nfilt);
    helimakelag(fac1,nx,nz);
    helimakelag(fac2,nx,nz);

    xkolmog_init(nkol);

    for (iw=0; iw < nw; iw++) { /* frequency loop */
	w = (iw+1)*dw;

	if (impresp) w=dw*nw/2;

	wov = w/v;
	wov2 = sf_cmplx(epsdamp,wov);
#ifdef SF_HAS_COMPLEX_H
	wov2 = -wov2*wov2;
#else
	wov2 = sf_cneg(sf_cmul(wov2,wov2));
#endif

	sf_warning("%g %g (%d of %d)",crealf(wov2),cimagf(wov2),iw,nw);

	init_wave(init,nx,dx,nz,dz,pp,wov,nw,iw);

	for (iz=0; iz < nx*nz; iz++) {
	    qq[iz]=sf_cmplx(0.,0.);
	}

	/* isotropic laplacian = 5-point laplacian */
	a= sf_cmplx(0.,0.);
#ifdef SF_HAS_COMPLEX_H
	b= gamma*(1+sixth*wov2)* (-1./(dz*dz));
	c= gamma*(1+sixth*wov2)* (-1./(dx*dx));
	d= gamma*(1+sixth*wov2)* (2/(dx*dx) + 2/(dz*dz))  -wov2;
#else
	b = sf_crmul(sf_cadd(sf_cmplx(1.,0.),sf_crmul(wov2,sixth)),
		     gamma*(-1./(dz*dz)));
	c = sf_crmul(sf_cadd(sf_cmplx(1.,0.),sf_crmul(wov2,sixth)),
		     gamma*(-1./(dx*dx)));
	d = sf_cadd(sf_crmul(sf_cadd(sf_cmplx(1.,0.),sf_crmul(wov2,sixth)),
			     gamma*(2/(dx*dx) + 2/(dz*dz))),sf_cneg(wov2));
#endif

	/* + rotated 5-point laplacian */
#ifdef SF_HAS_COMPLEX_H
	a += (1-gamma)*(1+sixth*wov2)* (-0.5/(dx*dz));
	b += (1-gamma)*(1+sixth*wov2)*0.;
	c += (1-gamma)*(1+sixth*wov2)*0.;
	d += (1-gamma)*(1+sixth*wov2)* 2.0/(dx*dz);
#else
	a = sf_cadd(a,sf_crmul(sf_cadd(sf_cmplx(1.0,0.0),
				       sf_crmul(wov2,sixth)),
			       (1-gamma)*(-0.5/(dx*dz))));
	d = sf_cadd(d,sf_crmul(sf_cadd(sf_cmplx(1.0,0.0),
				       sf_crmul(wov2,sixth)),
			       (1-gamma)*(2.0/(dx*dz))));
#endif

	aa->flt[0] = a; aa->lag[0] = -nx-1;
	aa->flt[1] = b; aa->lag[1] = -nx;
	aa->flt[2] = a; aa->lag[2] = -nx+1;

	aa->flt[3] = c; aa->lag[3] = -1;
	aa->flt[4] = d; aa->lag[4] = 0;
	aa->flt[5] = c; aa->lag[5] = 1;

	aa->flt[6] = a; aa->lag[6] = nx-1;
	aa->flt[7] = b; aa->lag[7] = nx;
	aa->flt[8] = a; aa->lag[8] = nx+1;

	xkolmog_helix(aa,fac1,fac2);

	for (i=0; i < nfilt; i++) {
#ifdef SF_HAS_COMPLEX_H
	    fac1->flt[i]=0.5*(fac2->flt[i]+conjf(fac1->flt[i]));
#else
	    fac1->flt[i]=sf_crmul(sf_cadd(fac2->flt[i],conjf(fac1->flt[i])),
				  0.5);
#endif
	}

	if (impresp) {
	    for (iz=0; iz < nx*nz; iz++) {
		pp[iz]=sf_cmplx(0.,0.);
	    }	    
	    pp[nx/2-1]=sf_cmplx(1.,0.);
	    sf_complexwrite(pp,2*nx,imp);
	}

	
	cpolydiv_init(nx*nz,fac2);
	cpolydiv_lop(false,false,nx*nz,nx*nz,pp,qq);

	if (impresp) {
	    sf_complexwrite(qq,2*nx,imp);
	    break;
	}

	/* back to time domain */
	for (it=0; it < nt; it++) {
	    cshift = cexpf(sf_cmplx( 0.,-w*it*dt));
	    for (ix=0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
		    ppp[it][ix][iz] += crealf(qq[ix+iz*nx]*cshift);
#else
		    ppp[it][ix][iz] += crealf(sf_cmul(qq[ix+iz*nx],cshift));
#endif
		}
	    }
	}

    } /* end frequency loop */

    sf_floatwrite(ppp[0][0],nt*nx*nz,out);

    exit(0);
}









