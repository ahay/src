/* 1-D lowrank wave propagation on staggered grid*/
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
#include "source.h"


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
    /*Flag*/
    bool verb, cmplx;        
    
    /*I/O*/
    sf_file Fsrc,Fo,Frec;    /* I/O files */
    sf_file left, right;
    sf_file Fvel, Fden, Ffft, Fic;
    
    /*I/O for MMS*/
    sf_file Fpsrc, Fvsrc, Fpint, Fvint;
    
    sf_axis at,ax;    /* cube axes */
    sf_axis icaxis;

    /* I/O arrays*/
    float *src; /*point source*/        
    float **lt, **rt;
    float *ic, *vel, *den, *c11;
    /*I/O arrays for MMS*/
    float **psrc, **vsrc, *pint, *vint;


    /*Grid index variables*/
    int it,im,ik,ix;    
    int nt, nx, m2, nk, nx2, n1, n2, pad1;
    float cx;
    float kx, dkx, kx0;
    float dx, dt, d1;
    
   
    sf_complex *cwavex, *cwavemx;
    float *record;
    float **wavex;
    float *curtxx, *pretxx;
    float *curvx, *prevx;
       
    /*source*/
    spara sp={0};
    int   srcrange;
    float srctrunc; 
    bool  srcmms, srcdecay;
    float slx=0;
    int   spx=0;
    
    /*options*/
    float gdep; int gp;
    bool  inject;
    int   it0, icnx;

    
    tstart = clock();
    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    
    /*model veloctiy & density*/
    Fvel = sf_input("vel");
    Fden = sf_input("den");

    /* setup I/O files */
    Fsrc = sf_input ("in");
    Fo = sf_output("out");
    Frec = sf_output ("rec");

    /*parameters of source*/
    if (!sf_getbool("srcmms", &srcmms)) srcmms = false;
    /*use MMS source */
    if (!srcmms && !sf_getfloat("slx", &slx)) sf_error("Need slx input");
    /*source location in x */
    if (!srcmms && slx<0.0) sf_error("slx need to be >=0.0");
    if (!sf_getbool("srcdecay", &srcdecay)) srcdecay=true;
    /*source decay*/
    if (!sf_getint("srcrange", &srcrange)) srcrange=10;
    /*source decay range*/
    if (!sf_getfloat("srctrunc", &srctrunc)) srctrunc=100;
    /*trunc source after srctrunc time (s)*/
    
    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if (!sf_getbool("inject", &inject)) inject=true; 
    /*=y inject source; =n initial condition*/
    if (!inject && srcmms) sf_error("Initial condition and MMS are conflicted");


    /* Read/Write axes */
    at = sf_iaxa(Fsrc,1); nt = sf_n(at); dt = sf_d(at);
    ax = sf_iaxa(Fvel,1); nx = sf_n(ax); dx = sf_d(ax);
     
    sf_oaxa(Fo,ax,1); 
    sf_oaxa(Fo,at,2);
    /*set for record*/
    sf_oaxa(Frec, at, 1);
    
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

    /*parameters of geometry*/
    if (!sf_getfloat("gdep", &gdep)) gdep = 0.0;
    /*depth of geophone (meter)*/
    if (gdep <0.0) sf_error("gdep need to be >=0.0");
    /*source and receiver location*/
    spx = (int)(slx/dx+0.5);
    gp  = (int)(gdep/dx+0.5);
    

    /* read source */
    src = sf_floatalloc(nt);
    sf_floatread(src,nt,Fsrc);
    
    /* Initial Condition */
    if (inject == false) {
	Fic     = sf_input("ic");  
	/*initial condition*/
	if (SF_FLOAT != sf_gettype(Fic)) sf_error("Need float input of ic");
	icaxis = sf_iaxa(Fic, 1); 
	icnx = sf_n(icaxis);
	if (nx != icnx) sf_error("I.C. and velocity should be the same size.");
	ic = sf_floatalloc(nx);
	sf_floatread(ic, nx, Fic);	
    } else {
	ic = NULL;
    }

    /* Method of Manufactured Solution*/
    if (inject && srcmms) {
	Fpsrc = sf_input("presrc");
	Fvsrc = sf_input("velsrc");
	Fpint = sf_input("preinit");
	Fvint = sf_input("velinit");
		
	if (SF_FLOAT != sf_gettype(Fpsrc)) sf_error("Need float input of presrc");
	if (SF_FLOAT != sf_gettype(Fvsrc)) sf_error("Need float input of velsrc");
	if (SF_FLOAT != sf_gettype(Fpint)) sf_error("Need float input of preinit");
	if (SF_FLOAT != sf_gettype(Fvint)) sf_error("Need float input of velinit");
	
	psrc = sf_floatalloc2(nx, nt);
	vsrc = sf_floatalloc2(nx, nt);
	pint = sf_floatalloc(nx);
	vint = sf_floatalloc(nx);
	
	sf_floatread(psrc[0], nx*nt, Fpsrc);
	sf_floatread(vsrc[0], nx*nt, Fvsrc);
	sf_floatread(pint, nx, Fpint);
	sf_floatread(vint, nx, Fvint);
    } else {
	psrc = NULL;
	vsrc = NULL;
	pint = NULL;
	vint = NULL;
    }
	
    curtxx  = sf_floatalloc(nx2);
    curvx   = sf_floatalloc(nx2);
    pretxx  = sf_floatalloc(nx);
    prevx   = sf_floatalloc(nx);
    

    cwavex  = sf_complexalloc(nk);
    cwavemx = sf_complexalloc(nk);
    wavex   = sf_floatalloc2(nx2,m2);
    
    record = sf_floatalloc(nt);
   
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
#ifdef SF_HAS_FFTW
	sf_warning("FFTW is defined");
#endif
#ifdef SF_HAS_COMPLEX_H
	sf_warning("Complex is defined");
#endif
	sf_warning("nx=%d dx=%f ", nx, dx);
	sf_warning("nk=%d dkx=%f kx0=%f", nk, dkx, kx0);
	sf_warning("nx2=%d", nx2);
	sf_warning("slx=%f, spx=%d, gdep=%f gp=%d",slx,spx,gdep,gp);
	sf_warning("======================================");
    } //End if

    /*set source*/
    sp.trunc=srctrunc;
    sp.srange=srcrange;
    sp.alpha=0.5;
    sp.decay=srcdecay?1:0; 

    /* MAIN LOOP */
    it0 = 0;
    if (inject == false) {
	it0 = 1;
	for(ix = 0; ix < nx; ix++) {
	    curtxx[ix] = ic[ix];
	    pretxx[ix] = curtxx[ix];
	}
	sf_floatwrite(curtxx,nx,Fo);
    }
    /* MMS */
    if (inject && srcmms ) {
	it0 = 0;
	for (ix=0; ix < nx; ix++) {
	    curtxx[ix] = pint[ix]; /*P(x,0)*/
	    pretxx[ix] = pint[ix];
	    prevx[ix]  = vint[ix]; /*U(x, -dt/2)*/
	}
    }
    
    for (it=it0; it<nt; it++) {
	if(verb) sf_warning("it=%d/%d;", it, nt-1);
	
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
	    curvx[ix] = -1*dt/den[ix]*cx + prevx[ix];  /*vx(t+dt/2) = -dt/rho*dp/dx(t) + vx(t-dt/2) */
	} //ix

	/* MMS */
	if (inject && srcmms) 
	    for (ix = 0; ix < nx; ix++)
		curvx[ix] += vsrc[it][ix]*dt;

	for (ix = 0; ix < nx; ix++)
	    prevx[ix] = curvx[ix];
	
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
	    
	    curtxx[ix] = -1*dt*c11[ix]*cx + pretxx[ix]; /*P(x,t+dt)*/
	}
	
	if (inject) {
	    if (!srcmms && (it*dt)<=sp.trunc ) {
		curtxx[spx] += src[it]*dt;
		//curtxx[spx] += 0.5*(srcp[it]+srcp[it+1])*dt;
		/*dp/dt(t+dt/2) = -rho*v^2*du/dx(t+dt/2) + s(t+dt)*/
	    }
	    if (srcmms)
		for (ix = 0; ix < nx; ix++) 
		    curtxx[ix] += psrc[it][ix]*dt;
	}
	/* write wavefield to output */
	sf_floatwrite(pretxx,nx,Fo);/* write P(x,t)*/
	record[it] = pretxx[gp];
	
	for (ix=0; ix<nx; ix++)
	    pretxx[ix] = curtxx[ix];
	 
    }/*End of MAIN LOOP*/
    
    if(verb) sf_warning(".");
    sf_floatwrite(record, nt, Frec);
    
    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsglr is: %f seconds << ", duration);
    exit (0);	
    
}


