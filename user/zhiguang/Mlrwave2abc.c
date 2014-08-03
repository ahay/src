/* Simple 2-D wave propagation */
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

#include "fft2.h"

int main(int argc, char* argv[])
{
    bool verb, cmplx;        
    int it,iz,im,ik,ix,i,j, snap;     /* index variables */
    int nt,nz,nx, m2, nk, nzx, nz2, nx2, nzx2, n2, pad1;
    float c, old, dt, z0, x0, dz, dx;
	int nb, rnz, rnx;
	float *bcx, *bcz, coef;

    float  *ww,*rr;      /* I/O arrays*/
    sf_complex *cwave, *cwavem;
    float **wave, *curr, *prev;

    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    float **lt, **rt;
    sf_file left, right, snaps;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
	if (!sf_getint("nb", &nb)) nb=40; /* padded boundary width */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fr,1); nz = sf_n(az); dz = sf_d(az); z0 = sf_o(az);
	rnz=nz-2*nb; sf_setn(az, rnz); sf_seto(az, z0+dz*nb);
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); dx = sf_d(ax); x0 = sf_o(ax);
	rnx=nx-2*nb; sf_setn(ax, rnx); sf_seto(ax, x0+dx*nb);

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    
    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */

    if (snap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	
	sf_oaxa(snaps,az,1); 
	sf_oaxa(snaps,ax,2);
	sf_oaxa(snaps,at,3);

	sf_putint(snaps,"n3",(nt-1)/snap+1);
	sf_putfloat(snaps,"d3",dt*snap);
	sf_putfloat(snaps,"o3",0.);
    } else {
	snaps = NULL;
    }

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
	if (!sf_getfloat("coef", &coef)) coef=0.01; /* decaying coefameter */

	/* Calculate the absorbing parameter */
	bcx=sf_floatalloc(nx);
	bcz=sf_floatalloc(nz);
	for(ix=0; ix<nx; ix++)
		bcx[ix]=1.0;
	for(iz=0; iz<nz; iz++){
		bcz[iz]=1.0;
	}
	for(ix=0; ix<nb; ix++){
		bcx[ix]=exp(-coef*coef*(nb-1-ix)*(nb-1-ix));
		bcz[ix]=bcx[ix];
		bcx[nx-1-ix]=bcx[ix];
		bcz[nz-1-ix]=bcz[ix];
	}

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

    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt);  sf_floatread(ww,nt ,Fw);
    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);

    curr = sf_floatalloc(nzx2);
    prev = sf_floatalloc(nzx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    ifft2_allocate(cwavem);

    for (iz=0; iz < nzx; iz++) {
	prev[iz]=0.;
    }

    for (iz=0; iz < nzx2; iz++) {
	curr[iz]=0.;
    }


    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* matrix multiplication */
	fft2(curr,cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_crmul(cwave[ik],rt[ik][im]);
#endif
	    }
	    ifft2(wave[im],cwavem);
	}

	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		
		old = c = curr[j];
		c += c + ww[it] * rr[i] - prev[i];
//		prev[i] = old;
        prev[i] = old*bcx[ix]*bcz[iz];

		for (im = 0; im < m2; im++) {
		    c += lt[im][i]*wave[im][j];
		}

		curr[j] = c*bcx[ix]*bcz[iz];
//   	curr[j] = c;
	    }
	    	
	    if (NULL != snaps && 0 == it%snap && ix>=nb && ix<nx-nb) {
		/* write wavefield snapshots */
		sf_floatwrite(curr+ix*nz2+nb,rnz,snaps);
	    }
	}
    }
    if(verb) sf_warning(".");   

    /* write final wavefield to output */
    for (ix = nb; ix < nx-nb; ix++) {
	sf_floatwrite(curr+ix*nz2+nb,rnz,Fo); 
    }
    
    exit (0);
}
