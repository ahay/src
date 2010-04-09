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
    bool verb;        
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx, nm,nk, nzx, nz2, nx2, nzx2, n2;
    float c, old;

    float  *ww,*rr;      /* I/O arrays*/
    sf_complex *cwave, *cwavem;
    float **wave, *curr, *prev;

    sf_file Fw,Fr,Fo;    /* I/O files */
    sf_axis at,az,ax;    /* cube axes */

    float **lt, **rt, **mm;
    sf_file left, right, mid;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); 
    az = sf_iaxa(Fr,1); nz = sf_n(az); 
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); 

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);
    
    nk = fft2_init(nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");
    mid = sf_input("mid");

    if (!sf_histint(mid,"n1",&nm)) sf_error("No n1= in mid");
    if (!sf_histint(mid,"n2",&n2) || n2 != nm) sf_error("Need n2=%d in mid",n2);

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&n2) || n2 != nm) sf_error("Need n2=%d in left",nm);
    
    if (!sf_histint(right,"n1",&n2) || n2 != nm) sf_error("Need n1=%d in right",nm);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
 
    lt = sf_floatalloc2(nzx,nm);
    mm = sf_floatalloc2(nm,nm);
    rt = sf_floatalloc2(nm,nk);

    sf_floatread(lt[0],nzx*nm,left);
    sf_floatread(rt[0],nm*nk,right);
    sf_floatread(mm[0],nm*nm,mid);

    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt);  sf_floatread(ww,nt ,Fw);
    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);

    curr = sf_floatalloc(nzx2);
    prev = sf_floatalloc(nzx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,nm);

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

	for (im = 0; im < nm; im++) {
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
		prev[i] = old;

		for (im = 0; im < nm; im++) {
		    for (ik = 0; ik < nm; ik++) {
			c += lt[ik][i]*mm[im][ik]*wave[im][j];
		    }
		}

		curr[j] = c;
	    }
	    	
	    /* write wavefield to output */
	    sf_floatwrite(curr+ix*nz2,nz,Fo);
	}
    }
    if(verb) sf_warning(".");    
    
    exit (0);
}
