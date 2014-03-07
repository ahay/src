/* Complex 2-D wave propagation (with kiss-fft)*/
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

#include "cfft2nsps.h"

int prop1(bool adj, sf_complex **ini, sf_complex **lt, sf_complex **rt, int nz, int nx, int nt, int m2, int nkzx, int pad1, sf_complex **cc)
{
    bool verb=true;
    /* index variables */
    int it,iz,ix,im,ik,i,j;
    int nz2,nx2,nk,nzx2;
    sf_complex c;
    /* wavefield */
    sf_complex **wave,**wave2, *curr, *currm, *cwave, *cwavem;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    nzx2 = nz2*nx2;

    if (nk!=nkzx) sf_error("nk discrepancy!");

    curr   = sf_complexalloc(nzx2);
    currm  = sf_complexalloc(nzx2);
    
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nk,m2);
    wave2  = sf_complexalloc2(nzx2,m2);

    if (adj) icfft2_allocate(cwave);
    else icfft2_allocate(cwavem);

    /* initialization */
    for (ix = 0; ix < nx2; ix++) {
	for (iz=0; iz < nz2; iz++) {
	    j = iz+ix*nz2;
	    if (ix<nx && iz<nz)
		curr[j] = ini[ix][iz];
	    else 
		curr[j] = sf_cmplx(0.,0.);
	}
    }

    /* MAIN LOOP */
    if (adj) { /*NSPS*/
	for (it=nt-1; it>-1; it--) {
	    if(verb) sf_warning("it=%d/%d;",it,nt);
	    
	    /* matrix multiplication */
	    for (im = 0; im < m2; im++) {
		for (ix = 0; ix < nx; ix++) {
		    for (iz=0; iz < nz; iz++) {
			i = iz+ix*nz;  /* original grid */
			j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
			currm[j] = conjf(lt[im][i])*curr[j];
#else
			currm[j] = sf_cmul(conjf(lt[im][i]), curr[j]);
#endif
		    }
		}
		cfft2(currm,wave[im]);
	    }
	    
	    for (ik = 0; ik < nk; ik++) {
		c = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += wave[im][ik]*conjf(rt[ik][im]);
#else
		    c += sf_cmul(wave[im][ik],conjf(rt[ik][im]));
#endif
		}
		cwave[ik] = c;
	    }
	    icfft2(curr,cwave);
	} /* time stepping */
    } else { 
	for (it=0; it<nt; it++) {
	    if(verb) sf_warning("it=%d/%d;",it,nt);
	    cfft2(curr,cwave);
	    
	    for (im = 0; im < m2; im++) {
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		    cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]);
#endif
		}
		icfft2(wave2[im],cwavem);
	    }
	    
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+ix*nz;  /* original grid */
		    j = iz+ix*nz2; /* padded grid */
		    c = sf_cmplx(0.,0.);
		    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
			c += lt[im][i]*wave2[im][j];
#else
			c += sf_cmul(lt[im][i], wave2[im][j]);
#endif
		    }
		    curr[j] = c;
		}
	    }
	} /* time stepping */
    }
    if(verb) sf_warning("."); 

    /* output final result*/
    for (ix=0; ix<nx; ix++)
	for (iz=0; iz<nz; iz++)
	    cc[ix][iz] = curr[iz+ix*nz2];
    
    cfft2_finalize();
    return 0;
}


int main(int argc, char* argv[])
{
    bool verb,adj;
    int nt,nz,nx,m2,nk,nzx,nz2,nx2,n2,pad1;
    float dt;

    sf_complex **ini, **cc;
    
    sf_file Fi,Fo;    /* I/O files */
    sf_axis az,ax;    /* cube axes */

    sf_complex **lt, **rt;
    sf_file left, right;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(!sf_getbool("adj",&adj)) adj=false; /* true: NSPS; default: PSPI */
    if(!sf_getint("nt",&nt)) sf_error("Need nt!"); /* verbosity */
    if(!sf_getfloat("dt",&dt)) sf_error("Need dt!"); /* verbosity */

    /* setup I/O files */
    Fi = sf_input ("in" );
    Fo = sf_output("out");

    sf_settype(Fo,SF_COMPLEX);

    /* Read/Write axes */
    az = sf_iaxa(Fi,1); nz = sf_n(az); 
    ax = sf_iaxa(Fi,2); nx = sf_n(ax); 

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 

    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

    nz2 = kiss_fft_next_fast_size(nz*pad1);
    nx2 = kiss_fft_next_fast_size(nx);
    nk = nz2*nx2; /*wavenumber*/
    nzx = nz*nx;

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
  
    lt = sf_complexalloc2(nzx,m2);
    rt = sf_complexalloc2(m2,nk);

    sf_complexread(lt[0],nzx*m2,left);
    sf_complexread(rt[0],m2*nk,right);

    sf_fileclose(left);
    sf_fileclose(right);

    ini=sf_complexalloc2(nz,nx);
    cc=sf_complexalloc2(nz,nx);

    sf_complexread(ini[0],nzx,Fi);
    sf_fileclose(Fi);

    /* wave propagation*/

    prop1(adj,ini, lt, rt, nz, nx, nt, m2, nk, 1, cc);

    /* output result */
    sf_complexwrite(cc[0], nzx, Fo);

    exit (0);
}

