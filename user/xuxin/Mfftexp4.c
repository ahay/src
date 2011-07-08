/* 2-D exploding reflector modeling/migration */

/*
  Copyright (C) 2011 KAUST
  
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
    bool mig,cmplx,verb;
    int it,ix,iz,nt,nx,nz,nx2,nz2,nzx,nzx2,nkx,nkz,n,im,i,j,m1,m2,ik,nk,pad1,jsnap;
    float d,c,old,*dd,**curr,**prev,**ll,**mm,**rr,**lm,**wave;
    sf_complex *cwave,*cwavem;
    sf_file Fi,Fo,Fu,Fl,Fm,Fr;
    sf_axis ax,az,at,aj;

    sf_init(argc,argv);

    if (!sf_getbool("mig",&mig)) mig=false;
    /* if n, modeling; if y, migration */
    if(!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getbool("cmplx",&cmplx)) cmplx=false;
    /* use complex FFT */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getint("pad1",&pad1)) pad1=1;
    /* padding factor on the first axis */ 
    if (!sf_getint("jsnap",&jsnap)) jsnap=1;
    /* jump in wavefield snapshots */
    if (!sf_getint("n",&n)) sf_error("Need n=");
    /* if mig=y, nz; if mig=n, nt */
    if (!sf_getfloat("d",&d)) sf_error("Need d=");
    /* if mig=y, dz; if mig=n, dt */

    Fi = sf_input("in");
    /* if mig=y, zero-offset setion; if mig=n, reflectivity */
    Fo = sf_output("out");
    /* if mig=y, migration image; if mig=n, zero-offset section */
    Fu = sf_output("wave");
    /* wavefield */

    /* axes */
    if (mig) {
	ax = sf_iaxa(Fi,1);
	at = sf_iaxa(Fi,2);
	az = sf_maxa(n,0,d);
	sf_setlabel(az,"z");
    } else {
	az = sf_iaxa(Fi,1);
	ax = sf_iaxa(Fi,2);
	at = sf_maxa(n,0,d);
	sf_setlabel(at,"t");
    }
    nx = sf_n(ax); 
    nz = sf_n(az); 
    nt = sf_n(at); 

    if (mig) {
	sf_oaxa(Fo,az,1);
	sf_oaxa(Fo,ax,2);
    } else {
	sf_oaxa(Fo,ax,1);
	sf_oaxa(Fo,at,2);
    }

    for (n=0, it=0; it < nt; it++)
	if (!(it % jsnap)) n++;
    aj = sf_maxa(n,sf_o(at),sf_d(at)*jsnap);
    sf_setlabel(aj,"t");
    sf_setunit(aj,"");
    sf_oaxa(Fu,az,1); 
    sf_oaxa(Fu,ax,2); 
    sf_oaxa(Fu,aj,3);

    /* wavenumbers */
    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);
    nzx = nz*nx;
    nzx2 = nz2*nx2;
    nkx = nx2;
    nkz = (int)(nk/nkx);
    if (verb) {
	sf_warning("nx =%d, nz =%d",nx,nz);
	sf_warning("nx2=%d, nz2=%d",nx2,nz2);
	sf_warning("nkx=%d, nkz=%d",nkx,nkz);
    }

    /* propagator matrices */
    Fl = sf_input("left");
    Fr = sf_input("right");
    Fm = sf_input("middle");

    if (!sf_histint(Fm,"n1",&m1)) sf_error("No n1= in middle");
    if (!sf_histint(Fm,"n2",&m2)) sf_error("No n2= in middle");
    if (verb)
	sf_warning("rank(k)=%d, rank(x)=%d",m1,m2);
    if (!sf_histint(Fl,"n1",&n) || n != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(Fl,"n2",&n) || n != m1)  sf_error("Need n2=%d in left",m1);
    if (!sf_histint(Fr,"n1",&n) || n != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(Fr,"n2",&n) || n != nk) sf_error("Need n2=%d in right",nk);

    ll = sf_floatalloc2(nzx,m1);
    mm = sf_floatalloc2(m1,m2);
    rr = sf_floatalloc2(m2,nk);
    
    sf_floatread(ll[0],nzx*m1,Fl);
    sf_floatread(mm[0],m1*m2,Fm);
    sf_floatread(rr[0],m2*nk,Fr);

    /* precompute left*middle */
    lm = sf_floatalloc2(nzx,m2);
    for (im=0; im < m2; im++) {
	for (ix=0; ix < nzx; ix++) {
	    for (c=0, ik=0; ik < m1; ik++)
		c += ll[ik][ix]*mm[im][ik];
	    lm[im][ix] = c;
	}
    }

    curr = sf_floatalloc2(nz2,nx2);
    prev = sf_floatalloc2(nz,nx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    for (ix=0; ix < nx*nz; ix++)
	prev[0][ix]=0.0f;

    for (ix=0; ix < nx2*nz2; ix++)
	curr[0][ix]=0.0f;

    /* read reflectivity */
    if (!mig) {
	for (ix=0; ix < nx; ix++) {
	    sf_floatread(curr[ix],nz,Fi);
	}
    }
    dd = sf_floatalloc(nx);

    /* MAIN LOOP */
    if (verb) sf_warning("nt=%d",nt);
    for (it=0; it < nt; it++) {
	sf_warning("it=%d;",(mig ? (nt-1-it) : it));

	/* inject source */
	if (mig) {
	    sf_floatread(dd,nx,Fi);
	    for (ix=0; ix < nx; ix++)
		curr[ix][0] += dd[ix];
	}
	
	/* matrix multiplication */
	fft2(curr[0],cwave);

	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rr[ik][im];
#else
		cwavem[ik] = sf_crmul(cwave[ik],rr[ik][im]);
#endif
	    }
	    ifft2(wave[im],cwavem);
	}

	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		
		old = curr[ix][iz];
		for (c = 0,im = 0; im < m2; im++) {
		    c += lm[im][i]*wave[im][j];
		}
		curr[ix][iz] = c - prev[ix][iz];
		prev[ix][iz] = old;
	    }

	    /* write wavefield */
	    if (!(it % jsnap))
		sf_floatwrite(curr[ix],nz,Fu);

	    /* extract data */
	    if (!mig)
		dd[ix] = curr[ix][0];
	}

	/* write data */
	if (!mig)
	    sf_floatwrite(dd,nx,Fo);	
    }
    if(verb) sf_warning(".");    

    /* write image */
    if (mig) {
	for (ix=0; ix < nx; ix++) {
	    /* for (iz=0; iz < nz; iz++)
		curr[ix][iz] *= -1; */
	    sf_floatwrite(curr[ix],nz,Fo);
	}
    }
    exit(0);
}
