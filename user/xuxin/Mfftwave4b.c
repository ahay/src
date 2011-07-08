/* attempting to achieve stable solution */
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
/* receivers are positioned at iz=0 grids.
   To put receivers at any coord, use lintr. 
   Not implemented atm */

#include <rsf.h>
#include <math.h>
#include "fft2.h"
#include "interp.h"

int main(int argc, char* argv[])
{
    bool verb,cmplx,wrtdat,stable;        
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx,m1,m2,nk,nzx,nz2,nx2,nzx2,n2,n,pad1,nkx,nkz,jsnap,ns;
    float c,old,*xs,*zs,*ww;
    lint2dp lints;

    sf_complex *cwave, *cwavem;
    float **wave,**curr,**prev, *dd;

    sf_file Fw,Fs,Fo,Fd,Fv;    /* I/O files */
    sf_axis at,az,ax,aj,px,pz;    /* cube axes */

    float **lt,**rt,**mm,**ltmm,*uu;
    sf_file left, right, mid;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(!sf_getbool("stable",&stable)) stable=false; /* stable correction */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fs = sf_input ("sou");
    Fv = sf_input ("vel"); /* read n,o,d only */
    Fo = sf_output("out");
    if (NULL != sf_getstring("dat")) {
	wrtdat = true;
	Fd = sf_output("dat");
    } else {
	wrtdat = false;
    }

    /* Read/Write axes */
    at = sf_iaxa(Fw,2); nt = sf_n(at); 
    az = sf_iaxa(Fv,1); nz = sf_n(az); 
    ax = sf_iaxa(Fv,2); nx = sf_n(ax); 

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if (!sf_getint("jsnap",&jsnap)) jsnap=1;
    /* jump in wavefield snapshots */

    for (n=0,it=0; it<nt; it++)
	if (!(it % jsnap)) n++;
    aj = sf_maxa(n,sf_o(at),sf_d(at)*jsnap);
    sf_setlabel(aj,"t");
    sf_setunit(aj,"");
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,aj,3);
    if (wrtdat) {
	sf_oaxa(Fd,ax,1);
	sf_oaxa(Fd,at,2);
	dd = sf_floatalloc(nx);
    }

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;
    nkx = nx2;
    nkz = (int)(nk/nkx);
    if (verb)
	sf_warning("nx=%d,nz=%d,nx=%d,nz2=%d,nkx=%d,nkz=%d",nx,nz,nx2,nz2,nkx,nkz);

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");
    mid = sf_input("middle");

    if (!sf_histint(mid,"n1",&m1)) sf_error("No n1= in mid");
    if (!sf_histint(mid,"n2",&m2)) sf_error("No n2= in mid");
    if (verb)
	sf_warning("rank(k)=%d, rank(x)=%d, nt=%d",m1,m2,nt);

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&n2) || n2 != m1)  sf_error("Need n2=%d in left",m1);
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lt = sf_floatalloc2(nzx,m1);
    mm = sf_floatalloc2(m1,m2);
    rt = sf_floatalloc2(m2,nk);
    
    sf_floatread(lt[0],nzx*m1,left);
    sf_floatread(mm[0],m1*m2,mid);
    sf_floatread(rt[0],m2*nk,right);

    /* stable trick */
    ltmm = sf_floatalloc2(nzx,m2); /* lt * mm */
    uu = sf_floatalloc(nzx); /* ltmm * rt[gd][..] */
    for (im=0; im < m2; im++) {
	for (ix=0; ix < nzx; ix++) {
	    c = 0;
	    for (ik=0; ik < m1; ik++)
		c += lt[ik][ix]*mm[im][ik];
	    ltmm[im][ix] = c;
	}
    }
    int gd=0; /* why pick k=0 ? */
    for (ix=0; ix < nzx; ix++) {
	c = 0;
	for (im=0; im < m2; im++) {
	    c += ltmm[im][ix]*rt[gd][im];
	}
	uu[ix] = abs(c) + SF_EPS;
    }
    
    if (stable) {
	for (ix=0; ix < nzx; ix++) {
	    for (im=0; im < m2; im++) {
		ltmm[im][ix] *= (uu[ix] > 2.0) ? 2.0/uu[ix] : 1.0;
		/* ltmm[im][ix] *= (uu[ix] > 0.01) ? 0.01/uu[ix] : 0.01; */
	    }
	}
    }

    /* source coord */
    if(!sf_histint(Fw,"n1",&ns)) sf_error("Need n1= in input");
    if (!sf_histint(Fs,"n1",&n) || n != ns) sf_error("Need n1=%d in sou",ns);
    xs = sf_floatalloc(ns);
    zs = sf_floatalloc(ns);
    ww = sf_floatalloc(ns);
    sf_floatread(xs,ns,Fs);
    sf_floatread(zs,ns,Fs);
    px = sf_maxa(nx2,sf_o(ax),sf_d(ax));
    pz = sf_maxa(nz2,sf_o(az),sf_d(az));
    lints = lint2d_init(zs,xs,ns,pz,px);
    curr = sf_floatalloc2(nz2,nx2);
    prev = sf_floatalloc2(nz,nx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    for (iz=0; iz < nzx; iz++) {
	prev[0][iz]=0.;
    }

    for (iz=0; iz < nzx2; iz++) {
	curr[0][iz]=0.;
    }

    /* MAIN LOOP */
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* inject source */
	sf_floatread(ww,ns,Fw);
	lint2d_inject(curr,ww,lints);

	/* matrix multiplication */
	fft2(curr[0],cwave);

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
		
		old = c = curr[ix][iz];
		c += c - prev[ix][iz];
		prev[ix][iz] = old;

		for (im = 0; im < m2; im++) {
		    c += ltmm[im][i]*wave[im][j];
		}

		curr[ix][iz] = c;
	    }
	    /* write wavefield to output */
	    if (!(it % jsnap))
		sf_floatwrite(curr[ix],nz,Fo);
 	}

	/* write data */
	if (wrtdat) {
	    for (ix=0; ix<nx; ix++)
		dd[ix] = curr[ix][0];
	    sf_floatwrite(dd,nx,Fd);
	}
    }
    if(verb) sf_warning(".");    
    
    exit (0);
}
