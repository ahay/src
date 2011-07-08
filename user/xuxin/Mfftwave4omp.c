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
/* revision 185 is working. This one attempts to improve speed by omp */
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "fft2.h"

int main(int argc, char* argv[])
{
    bool verb, cmplx;        
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nt,nz,nx,m1,m2,nk,nzx,nz2,nx2,nzx2,n2,n,pad1,nkx,nkz,jsnap,ompnth,ompath,ompchunk;
    float c, old;

    float  *ww,*rr;      /* I/O arrays*/
    sf_complex *cwave, *cwavem;
    float **wave,*curr,*prev,*dd,*lmr;

    sf_file Fw,Fr,Fo,Fd;    /* I/O files */
    sf_axis at,az,ax,aj;    /* cube axes */

    float **lt, **rt, **mm;
    sf_file left, right, mid;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fr = sf_input ("ref");
    Fo = sf_output("out");
    Fd = sf_output("dat");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); 
    az = sf_iaxa(Fr,1); nz = sf_n(az); 
    ax = sf_iaxa(Fr,2); nx = sf_n(ax); 

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
    if (!sf_getint("jsnap",&jsnap)) jsnap=1;
    /* jump in wavefield snapshots */
    if (!sf_getint("ompnth",&ompnth)) ompnth=1;
    /* OpenMP num of threads */
    if (!sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    /* OpenMP chunksize */
#ifdef _OPENMP
#pragma omp parallel
    ompath=omp_get_num_threads();
    omp_set_num_threads(ompnth);
    if (verb)
	sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    for (n=0,it=0; it<nt; it++)
	if (!(it % jsnap)) n++;
    aj = sf_maxa(n,sf_o(at),sf_d(at)*jsnap);
    sf_setlabel(aj,"t");
    sf_setunit(aj,"");
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,aj,3);
    sf_oaxa(Fd,ax,1);
    sf_oaxa(Fd,at,2);

    nk = fft2_init(cmplx,pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;
    nkx = nx2;
    nkz = (int)(nk/nkx);
    if (verb)
	sf_warning("nx =%4d,nz =%4d,nx2=%4d,nz2=%4d,nkx=%4d,nkz=%4d",nx,nz,nx2,nz2,nkx,nkz);

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");
    mid = sf_input("middle");

    if (!sf_histint(mid,"n1",&m1)) sf_error("No n1= in mid");
    if (!sf_histint(mid,"n2",&m2)) sf_error("No n2= in mid");
    if (verb)
	sf_warning("rank(k)=%d, rank(x)=%d",m1,m2);

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

    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt);  sf_floatread(ww,nt ,Fw);
    rr=sf_floatalloc(nzx); sf_floatread(rr,nzx,Fr);

    curr = sf_floatalloc(nzx2);
    prev = sf_floatalloc(nzx);
    dd   = sf_floatalloc(nx);
    lmr  = sf_floatalloc(nzx);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nzx2,m2);

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

#ifdef OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,i,j,im,ik) \
    share(nx,nz,nz2,m2,m1,lmr,lt,mm,wave)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		for (im = 0; im < m2; im++) {
		    for (ik = 0; ik < m1; ik++) {
			lmr[i] = lt[ik][i]*mm[im][ik]*wave[im][j];
		    }
		}
	    }
	}
#ifdef OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ix,iz,i,j,old,c)			    \
    share(nx,nz,nz2,prev,curr,ww,it,rr,lmr)
#endif
	for (ix = 0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		i = iz+ix*nz;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		
		old = c = curr[j];
		c += c + ww[it] * rr[i] - prev[i] + lmr[i];
		prev[i] = old;
		curr[j] = c;
	    }
	    /* write wavefield to output */
	    if (!(it % jsnap))
		sf_floatwrite(curr+ix*nz2,nz,Fo);
 	}

	/* write data */
	for (ix=0; ix<nx; ix++)
	    dd[ix] = curr[ix*nz2];
	sf_floatwrite(dd,nx,Fd);
    }
    if(verb) sf_warning(".");    
    
    exit (0);
}
