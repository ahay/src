/* Simple 2-D wave propagation */

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
#include "interp.h"

int main(int argc, char* argv[])
{
    bool verb,cmplx;
    int it,iz,im,ik,ix,i,j,nt,nz,nx,m1,m2,nk,nzx,nz2,nx2,nzx2,n,pad1,nkx,nkz,jsnap,ns,nr;
    float c,old,ox,dx,oz,dz,*xs,*zs,*xr,*zr,*ww;
    lint2dp lints,lintr;

    sf_complex *cwave,*cwavem;
    float **wave,**curr,**prev,*dd;

    sf_file Fi,Fss,Frr,Fo,Fu,Fv;
    sf_axis at,az,ax,aj,px,pz;

    float **ll,**rr,**mm,**lm;
    sf_file Fl,Fr,Fm;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getbool("cmplx",&cmplx)) cmplx=false;
    /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1;
    /* padding factor on the first axis */
    if (!sf_getint("jsnap",&jsnap)) jsnap=1;
    /* jump in wavefield snapshots */

    /* files */
    Fi = sf_input ("in" );
    Fss= sf_input ("sou");
    Fv = sf_input ("vel"); /* read n,o,d only */
    Fo = sf_output("out");
    Fu = sf_output("wave");

    /* axes */
    at = sf_iaxa(Fi,2); nt = sf_n(at); 
    az = sf_iaxa(Fv,1); nz = sf_n(az);
    ax = sf_iaxa(Fv,2); nx = sf_n(ax);
    ox = sf_o(ax); dx = sf_d(ax);
    oz = sf_o(az); dz = sf_d(az);

    for (n=0,it=0; it<nt; it++)
	if (!(it % jsnap)) n++;
    aj = sf_maxa(n,sf_o(at),sf_d(at)*jsnap);
    sf_setlabel(aj,"t");
    sf_setunit(aj,"");
    sf_oaxa(Fu,az,1); 
    sf_oaxa(Fu,ax,2); 
    sf_oaxa(Fu,aj,3);

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

    /* source coord */
    if (!sf_histint(Fi,"n1",&ns)) sf_error("Need n1= in input");
    if (!sf_histint(Fss,"n1",&n) || n != ns) sf_error("Need n1=%d in sou",ns);
    xs = sf_floatalloc(ns);
    zs = sf_floatalloc(ns);
    ww = sf_floatalloc(ns);
    sf_floatread(xs,ns,Fss);
    sf_floatread(zs,ns,Fss);
    px = sf_maxa(nx2,ox,dx);
    pz = sf_maxa(nz2,oz,dz);
    lints = lint2d_init(zs,xs,ns,pz,px);

    /* receiver coord */
    if (NULL != sf_getstring("rec")) {
	Frr = sf_input("rec");
	if (!sf_histint(Frr,"n1",&nr)) sf_error("Need n1= in rec");
	xr = sf_floatalloc(nr);
	zr = sf_floatalloc(nr);
	sf_floatread(xr,nr,Frr);
	sf_floatread(zr,nr,Frr);
    } else {
	nr = nx;
	xr = sf_floatalloc(nr);
	zr = sf_floatalloc(nr);
	for (xr[0]=ox,ix=1; ix < nr; ix++) {
	    xr[ix] = xr[ix-1] + dx;
	    zr[ix] = 0;
	}
    }
    lintr = lint2d_init(zr,xr,nr,pz,px);
    sf_putint(Fo,"n1",nr);
    sf_oaxa(Fo,at,2);

    curr = sf_floatalloc2(nz2,nx2);
    prev = sf_floatalloc2(nz,nx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

    for (iz=0; iz < nzx; iz++)
	prev[0][iz]=0.0f;   

    for (iz=0; iz < nzx2; iz++)
	curr[0][iz]=0.0f;   

    dd = sf_floatalloc(nr);

    /* MAIN LOOP */
    if (verb) sf_warning("nt=%d",nt);
    for (it=0; it<nt; it++) {
	if(verb) sf_warning("it=%d;",it);

	/* inject source */
	sf_floatread(ww,ns,Fi);
	lint2d_inject(curr,ww,lints);

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
	    lint2d_extract(curr,dd,lintr);
 	}

	/* write data */
	sf_floatwrite(dd,nr,Fo);
    }
    if(verb) sf_warning(".");    
    
    exit (0);
}
