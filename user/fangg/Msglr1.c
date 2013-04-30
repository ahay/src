/* Simple 1-D wave propagation on staggered grid*/
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
#include <time.h>

#include "fft1.h"
//#include "fft2.h"

sf_complex fplus(float _kx, float dx)
/*i*kx*exp(i*kx*dx/2)*/
{
    sf_complex res;
    float kx=_kx*2*SF_PI;
    float r = -kx*sinf(kx*dx*0.5);
    float i = kx*cosf(kx*dx*0.5);
    res = sf_cmplx(r,i);
    return res;
}

sf_complex fminu(float _kx, float dx)
/*i*kx*exp(-i*kx*dx/2)*/
{
    sf_complex res;
    float kx=_kx*2*SF_PI;
    float i = kx*cosf(kx*dx*0.5);
    float r = kx*sinf(kx*dx*0.5);
    res = sf_cmplx(r,i);
    return res;
}

int main(int argc, char* argv[])
{
    clock_t tstart, tend;
    double duration;
    bool verb, cmplx;        
    int it,im,ik,ix;     /* index variables */
    int nt, nx, m2, nk, nx2, n1, n2, pad1;
    float cx;
    float kx, dkx, kx0;
    float dx, dt, d1;
    
    float  *ww,*rr;      /* I/O arrays*/
    sf_complex *cwavex, *cwavemx;
    
    float **wavex;
    float *curtxx, *pretxx;
    float *curvx, *prevx;
    

    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,ax;    /* cube axes */

    float **lt, **rt;
    sf_file left, right;

    float *vel, *den, *c11;
    sf_file Fvel, Fden, Ffft, Fic;
    
    bool inject;
    float *ic;
    int icnx;
    sf_axis icaxis;
    
    tstart = clock();
    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */

    /* setup I/O files */
    Fw = sf_input ("in");
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    ax = sf_iaxa(Fr,1); nx = sf_n(ax); dx = sf_d(ax);
    
    
    sf_oaxa(Fo,ax,1); 
    sf_oaxa(Fo,at,2);

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if (!sf_getbool("inject", &inject)) inject=true; 
    /*=y inject source; =n initial condition*/

    nk = fft1_init(nx,&nx2);

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nx) sf_error("Need n1=%d in left",nx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");

    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);

    lt = sf_floatalloc2(nx,m2);
    rt = sf_floatalloc2(m2,nk);

    sf_floatread(lt[0],nx*m2,left);
    sf_floatread(rt[0],m2*nk,right);

    /*model veloctiy & density*/
    Fvel = sf_input("vel");
    Fden = sf_input("den");
    
    if (!sf_histint(Fvel,"n1", &n1) || n1 != nx) sf_error("Need n1=%d in vel", nx);
    if (!sf_histfloat(Fvel,"d1", &d1) || d1 != dx) sf_error("Need d1=%d in vel", dx);
    
    if (!sf_histint(Fden,"n1", &n1) || n1 != nx) sf_error("Need n1=%d in den", nx);
    if (!sf_histfloat(Fden,"d1", &d1) || d1 != dx) sf_error("Need d1=%d in den", dx);
    
    vel = sf_floatalloc(nx);
    den = sf_floatalloc(nx);
    c11 = sf_floatalloc(nx);

    sf_floatread(vel, nx, Fvel);
    sf_floatread(den, nx, Fden);
    
    for (ix = 0; ix < nx; ix++) {
	c11[ix] = den[ix]*vel[ix]*vel[ix];
	if (c11[ix] == 0) sf_warning("C11[%d] = %f",ix, c11[ix]);
    }
    
    /*parameters of fft*/
    Ffft = sf_input("fft");
    if (!sf_histint(Ffft,"n1", &nk)) sf_error("Need n1 in fft");
    if (!sf_histfloat(Ffft,"d1", &dkx)) sf_error("Need d1 in fft");
    if (!sf_histfloat(Ffft,"o1", &kx0)) sf_error("Need o1 in fft");
    
    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt); sf_floatread(ww, nt, Fw);
    rr=sf_floatalloc(nx); sf_floatread(rr,nx,Fr);

    curtxx = sf_floatalloc(nx2);
    curvx  = sf_floatalloc(nx2);
    pretxx  = sf_floatalloc(nx);
    prevx   = sf_floatalloc(nx);
    

    cwavex = sf_complexalloc(nk);
    cwavemx = sf_complexalloc(nk);
    wavex = sf_floatalloc2(nx2,m2);

    /*Initial Condition*/
    if (inject == false) {
	Fic     = sf_input("ic");  
	/*initial condition*/
	if (SF_FLOAT != sf_gettype(Fic)) sf_error("Need float input of ic");
	icaxis = sf_iaxa(Fic, 1); 
	icnx = sf_n(icaxis);
	if (nx != icnx) sf_error("I.C. and velocity should be the same size.");
	ic = sf_floatalloc(nx);
	sf_floatread(ic, nx, Fic);	
    }

    
    ifft1_allocate(cwavemx);
    
    for (ix=0; ix < nx; ix++) {
	pretxx[ix]=0.;
	prevx[ix] =0.;
    }

    for (ix=0; ix < nx2; ix++) {
	curtxx[ix]=0.;
	curvx[ix]=0.;
    }

    /* Check parameters*/
    if(verb) {
	sf_warning("======================================");
	sf_warning("nx=%d dx=%f ", nx, dx);
	sf_warning("nk=%d dkx=%f kx0=%f", nk, dkx, kx0);
	sf_warning("nx2=%d", nx2);
	sf_warning("======================================");
    } //End if
   
    /* MAIN LOOP */
   for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;", it);
	
	/*vx--- matrix multiplication */
	fft1(curtxx,cwavex);   /* P^(k,t) */
	
	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
		kx = kx0+dkx*ik;
#ifdef SF_HAS_COMPLEX_H
		cwavemx[ik] = cwavex[ik]*rt[ik][im];
		cwavemx[ik] = fplus(kx,dx)*cwavemx[ik];
#else
		cwavemx[ik] = sf_crmul(cwavex[ik],rt[ik][im]);
		cwavemx[ik] = sf_cmul(fplus(kx,dx), cwavex[ik]);
#endif
		
	    } // ik
	    ifft1(wavex[im], cwavemx); /* dp/dx  */
	} // im
	
	for (ix = 0; ix < nx; ix++) {
	    cx = 0.0;
	    for (im=0; im<m2; im++) {
		cx += lt[im][ix]*wavex[im][ix];
	    } 
	    curvx[ix] = -1*dt/den[ix]*cx + prevx[ix];  /*vx(t+dt) = -dt/rho*dp/dx + vx(t) */
	    prevx[ix] = curvx[ix];
	} //ix
	
	/*txx--- matrix multiplication */
	fft1(curvx, cwavex);
	
	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++ ) {
		kx = kx0 + dkx*ik;
		
#ifdef SF_HAS_COMPLEX_H
		cwavemx[ik] = cwavex[ik]*rt[ik][im];
		cwavemx[ik] = fminu(kx,dx)*cwavemx[ik];
#else
		cwavemx[ik] = sf_crmul(cwavex[ik],rt[ik][im]);
		cwavemx[ik] = sf_cmul(fminu(kx,dx), cwavemx[ik]);
		
#endif
	    }
	    ifft1(wavex[im], cwavemx); /* dux/dx */
	}
	
	for (ix = 0; ix < nx; ix++) {
	    cx = 0.0;
	    for (im=0; im<m2; im++) {
		cx += lt[im][ix]*wavex[im][ix];
	    }
	    if (inject == true) {
		curtxx[ix] = -1*dt*c11[ix]*cx + pretxx[ix]+ ww[it] * rr[ix];
	    } else {
		curtxx[ix] = -1*dt*c11[ix]*cx + pretxx[ix];
	    }

	    pretxx[ix] = curtxx[ix];
	    
	    /* write wavefield to output */
	    
	}

	 
	if (it==0 && inject == false) {
	    for(ix = 0; ix < nx; ix++) {
		curtxx[ix] = ic[ix];
		pretxx[ix] = curtxx[ix];
		//vxn0[ix+marg+pmlout] = ic[ix];
	    }
	}
		
	sf_floatwrite(curtxx,nx,Fo);
	
   }/*End of MAIN LOOP*/
   
   if(verb) sf_warning(".");
   
   tend = clock();
   duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
   sf_warning(">> The CPU time of sfsglr is: %f seconds << ", duration);
   exit (0);	
   
}


