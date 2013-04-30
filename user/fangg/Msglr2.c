/* Simple 2-D wave propagation on staggered grid*/
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

#include "fft2.h"

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
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx, m2, nk, nkx, nkz, nzx, nz2, nx2, nzx2, n1, n2, pad1;
    float cx, cz;
    float kx, kz, dkx, dkz, kx0, kz0;
    float dx, dz, dt, d1, d2;
    
    float  *ww,*rr;      /* I/O arrays*/
    sf_complex *cwavex, *cwavez, *cwavemx, *cwavemz;
    float **wavex, **wavez;
    float *curtxx, *pretxx;
    float *curvx, *prevx, *curvz, *prevz;
    

    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    float **lt, **rt;
    sf_file left, right;

    float **vel, **den, **c11;
    sf_file Fvel, Fden, Ffft; 

    tstart = clock();
    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fr,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); dx = sf_d(ax);
    
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");

    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);

    lt = sf_floatalloc2(nzx,m2);
    rt = sf_floatalloc2(m2,nk);

    sf_floatread(lt[0],nzx*m2,left);
    sf_floatread(rt[0],m2*nk,right);

    /*model veloctiy & density*/
    Fvel = sf_input("vel");
    Fden = sf_input("den");
    
    if (!sf_histint(Fvel,"n1", &n1) || n1 != nz) sf_error("Need n1=%d in vel", nz);
    if (!sf_histfloat(Fvel,"d1", &d1) || d1 != dz) sf_error("Need d1=%d in vel", dz);
    if (!sf_histint(Fvel,"n2", &n2) || n2 != nx) sf_error("Need n2=%d in vel", nx);
    if (!sf_histfloat(Fvel,"d2", &d2) || d2 != dx) sf_error("Need d2=%d in vel", dx);

    if (!sf_histint(Fden,"n1", &n1) || n1 != nz) sf_error("Need n1=%d in den", nz);
    if (!sf_histfloat(Fden,"d1", &d1) || d1 != dz) sf_error("Need d1=%d in den", dz);
    if (!sf_histint(Fden,"n2", &n2) || n2 != nx) sf_error("Need n2=%d in den", nx);
    if (!sf_histfloat(Fden,"d2", &d2) || d2 != dx) sf_error("Need d2=%d in den", dx);

    vel = sf_floatalloc2(nz, nx);
    den = sf_floatalloc2(nz, nx);
    c11 = sf_floatalloc2(nz, nx);

    sf_floatread(vel[0], nzx, Fvel);
    sf_floatread(den[0], nzx, Fden);
    
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    c11[ix][iz] = den[ix][iz]*vel[ix][iz]*vel[ix][iz];
	}
    }
    
    /*parameters of fft*/
    Ffft = sf_input("fft");
    if (!sf_histint(Ffft,"n1", &nkz)) sf_error("Need n1 in fft");
    if (!sf_histint(Ffft,"n2", &nkx)) sf_error("Need n2 in fft");
    if ( nkx*nkz != nk )  sf_error("Need nk=nkx*nkz, nk=%d, nkx=%d, nkz=%d", nk, nkx, nkz);
    if (!sf_histfloat(Ffft,"d1", &dkz)) sf_error("Need d1 in fft");
    if (!sf_histfloat(Ffft,"d2", &dkx)) sf_error("Need d2 in fft");
    if (!sf_histfloat(Ffft,"o1", &kz0)) sf_error("Need o1 in fft");
    if (!sf_histfloat(Ffft,"o2", &kx0)) sf_error("Need o2 in fft");

    

    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt); sf_floatread(ww, nt, Fw);
    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);

    curtxx = sf_floatalloc(nzx2);
    curvx  = sf_floatalloc(nzx2);
    curvz  = sf_floatalloc(nzx2);
    pretxx  = sf_floatalloc(nzx);
    prevx   = sf_floatalloc(nzx);
    prevz   = sf_floatalloc(nzx);

    cwavex = sf_complexalloc(nk);
    cwavez = sf_complexalloc(nk);
    cwavemx = sf_complexalloc(nk);
    cwavemz = sf_complexalloc(nk);
    wavex = sf_floatalloc2(nzx2,m2);
    wavez = sf_floatalloc2(nzx2,m2);

    ifft2_allocate(cwavemx);
    ifft2_allocate(cwavemz);

    for (iz=0; iz < nzx; iz++) {
	pretxx[iz]=0.;
	prevx[iz] =0.;
	prevz[iz] =0.;
    }

    for (iz=0; iz < nzx2; iz++) {
	curtxx[iz]=0.;
	curvx[iz]=0.;
	curvz[iz]=0.;
    }

    /* Check parameters*/
    if(verb) {
	sf_warning("======================================");
	sf_warning("nx=%d nz=%d nzx=%d dx=%f dz=%f", nx, nz, nzx, dx, dz);
	sf_warning("nkx=%d nkz=%d dkx=%f dkz=%f nk=%d", nkx, nkz, dkx, dkz, nk);
	sf_warning("nx2=%d nz2=%d nzx2=%d", nx2, nz2, nzx2);
	sf_warning("======================================");
    } //End if
    
   
    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);
	
	/*vx, vz--- matrix multiplication */
	fft2(curtxx,cwavex);   /* P^(k,t) */
	
	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
		kx = kx0+dkx*(ik/nkz);  
		kz = kz0+dkz*(ik%nkz);
		
#ifdef SF_HAS_COMPLEX_H
		cwavemz[ik] = cwavex[ik]*rt[ik][im];
		cwavemx[ik] = fplus(kx,dx)*cwavemz[ik];
		cwavemz[ik] = fplus(kz,dz)*cwavemz[ik];
#else
		cwavemz[ik] = sf_crmul(cwavex[ik],rt[ik][im]);
		cwavemx[ik] = sf_cmul(fplus(kx,dx), cwavemz[ik]);
		cwavemz[ik] = sf_cmul(fplus(kz,dz), cwavemz[ik]);
#endif
	    }
	    ifft2(wavex[im], cwavemx); /* dp/dx  */
	    ifft2(wavez[im], cwavemz); /* dp/dz  */
	    
	}
	
	for (ix = 0; ix < nx; ix++) {
	    for (iz = 0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */

		cx = 0.0;
		cz = 0.0;
		for (im=0; im<m2; im++) {
		    cx += lt[im][i]*wavex[im][j];
		    cz += lt[im][i]*wavez[im][j];
		}

		curvx[j] = -1*dt/den[ix][iz]*cx + prevx[i];  /*vx(t+dt) = -dt/rho*dp/dx + vx(t) */
		prevx[i] = curvx[j];

		curvz[j] = -1*dt/den[ix][iz]*cz + prevz[i];
		prevz[i] = curvz[j];
	    }
	}
	
	/*txx--- matrix multiplication */
	fft2(curvx, cwavex);
	fft2(curvz, cwavez);
	
	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++ ) {
		kx = kx0 + dkx*(ik/nkz);
		kz = kz0 + dkz*(ik%nkz);

#ifdef SF_HAS_COMPLEX_H
		cwavemz[ik] = cwavez[ik]*rt[ik][im];
		cwavemx[ik] = cwavex[ik]*rt[ik][im];
		cwavemx[ik] = fminu(kx,dx)*cwavemx[ik];
		cwavemz[ik] = fminu(kz,dz)*cwavemz[ik];
#else
		cwavemz[ik] = sf_crmul(cwavez[ik],rt[ik][im]);
		cwavemx[ik] = sf_crmul(cwavex[ik],rt[ik][im]);
		cwavemx[ik] = sf_cmul(fplus(kx,dx), cwavemx[ik]);
		cwavemz[ik] = sf_cmul(fplus(kz,dz), cwavemz[ik]);
#endif
	    }
	    ifft2(wavex[im], cwavemx); /* dux/dx */
	    ifft2(wavez[im], cwavemz); /* duz/dz */
	}

	for (ix = 0; ix < nx; ix++) {
	    for (iz = 0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		
		cx = 0.0;
		cz = 0.0;
		
		for (im=0; im<m2; im++) {
		    cx += lt[im][i]*wavex[im][j];
		    cz += lt[im][i]*wavez[im][j];
		}	
		
		curtxx[j] = -1*dt*c11[ix][iz]*(cx+cz) + pretxx[i]+ ww[it] * rr[i];
		pretxx[i] = curtxx[j];
	    }
	    
	    /* write wavefield to output */
	    sf_floatwrite(curtxx+ix*nz2,nz,Fo); 
	}

    }/*End of MAIN LOOP*/

    if(verb) sf_warning(".");
    
    tend = clock();
    duration=(double)(tend-tstart)/CLOCKS_PER_SEC;
    sf_warning(">> The CPU time of sfsglr is: %f seconds << ", duration);
    exit (0);	

}






    
 

