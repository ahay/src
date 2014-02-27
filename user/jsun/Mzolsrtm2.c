/* 2-D FFT-based zero-offset exploding reflector modeling/migration  */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cfft2nsps.h"
#include "timer.h"
/*******************************************************/
/* wave propagation utils*/
typedef struct Geopar {
    /*geometry parameters*/
    int   nx;
    int   nz;
    float dx;
    float dz;
    float ox;
    int   gpz;
    /*time parameters*/
    int nt;
    float dt;
    int snap;
    int wfnt;
    /*prop matrix*/
    int nzx2;
    int nk;
    int m2;
    /*fft parameters*/
    int nkz;
    int nkx;
    float dkz;
    float dkx;
    float kz0;
    float kx0;
    /*lowrank parameters*/
    int m; 
    int n2; // left: m*n2
    int n;  //right: n2*n
} * geopar; /*geometry parameters*/

int lrexp(sf_complex **img, sf_complex **dat, bool adj, sf_complex **lt, sf_complex **rt, geopar geop, int pad1, bool verb, int snap, sf_complex ***wvfld)
/*< zero-offset exploding reflector modeling/migration >*/
{
    int it, nt, ix, nx, nx2, iz, nz, nz2, nzx2, gpz, wfnt, wfit;
    int im, i, j, m2, ik, nk;
    float dt, dx, dz, ox;
    sf_complex *curr, **wave, *cwave, *cwavem, c;
    sf_complex *currm;

    nx  = geop->nx;
    nz  = geop->nz;
    dx  = geop->dx;
    dz  = geop->dz;
    ox  = geop->ox;
    gpz = geop->gpz;
    nt  = geop->nt;
    dt  = geop->dt;
    snap= geop->snap;
    nzx2= geop->nzx2;
    m2  = geop->m2;
    wfnt= geop->wfnt;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    if (nk!=geop->nk) sf_error("nk discrepancy!");

    curr = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);
    if (adj) {
	currm  = sf_complexalloc(nzx2);
	icfft2_allocate(cwave);
    } else {
	cwavem = sf_complexalloc(nk);
	icfft2_allocate(cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

    if (adj) { /* migration <- read data */
	wfit = (int)(nt-1)/snap; // wfnt-1
	/* time stepping */
	for (it=nt-1; it > -1; it--) {
	    if (verb) sf_warning("it=%d;",it);
	
	    /* matrix multiplication */
	    for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j) shared(currm,lt,curr)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for private(ik,im,c) shared(wave,rt,cwave)
#endif	    
	    for (ik = 0; ik < nk; ik++) {
		c = sf_cmplx(0.,0.);
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += wave[im][ik]*conjf(rt[ik][im]);
#else
		    c += sf_cmul(wave[im][ik],conjf(rt[ik][im])); //complex multiplies complex
#endif
		}
		cwave[ik] = c;
	    }

	    icfft2(curr,cwave);

#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    /* data injection */
	    for (ix=0; ix < nx; ix++) {
		curr[gpz+ix*nz2] += dat[ix][it];
	    }

	    if (snap > 0 && it%snap == 0) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
		for ( ix = 0; ix < nx; ix++) {
		    for ( iz = 0; iz<nz; iz++ ) { 
			j = iz+ix*nz2; /* padded grid */
			wvfld[wfit][ix][iz] = curr[j];
		    }
		}
		wfit--;
	    }
	} /*time iteration*/
	/*generate image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		img[ix][iz] = curr[iz+ix*nz2];
	    }
	}
    } else { /* modeling -> write data */
	/*exploding reflector*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		curr[iz+ix*nz2]=img[ix][iz];
	    }
	}
	wfit = 0;
	/* time stepping */
	for (it=0; it < nt; it++) {
	    if (verb) sf_warning("it=%d;",it);
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    /* record data */
	    for (ix=0; ix < nx; ix++) {
		dat[ix][it] = curr[gpz+ix*nz2];
	    }
	    /* matrix multiplication */
	    cfft2(curr,cwave);
	    
	    for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		    cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]);
#endif
		}
		icfft2(wave[im],cwavem);
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j,im,c) shared(curr,lt,wave)
#endif
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+ix*nz;  /* original grid */
		    j = iz+ix*nz2; /* padded grid */
		    
		    c = sf_cmplx(0.,0.); /* initialize */
		    
		    for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
			c += lt[im][i]*wave[im][j];
#else
			c += sf_cmul(lt[im][i], wave[im][j]);
#endif	    
		    }
		    curr[j] = c;
		}
	    }
	    if (snap > 0 && it%snap == 0) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
		for ( ix = 0; ix < nx; ix++) {
		    for ( iz = 0; iz<nz; iz++ ) { 
			j = iz+ix*nz2; /* padded grid */
			wvfld[wfit][ix][iz] = curr[j];
		    }
		}
		wfit++;
	    }
	}
    }
    if (verb) sf_warning(".");

    return 0;
}

int main(int argc, char* argv[])
{
    bool adj,timer,verb;
    int nt, nx, nz, nx2, nz2, nzx, nzx2, ntx, pad1, snap, gpz, wfnt;
    int m2, n2, nk, nth=1;
    float dt, dx, dz, ox;
    sf_complex **img, **dat, **lt, **rt, ***wvfld;
    sf_file data, image, left, right, snaps;
    double time=0.,t0=0.,t1=0.;
    geopar geop;

    sf_init(argc,argv);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if n, modeling; if y, migration */
    if(! sf_getbool("timer",&timer)) timer=false;
    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */
    if (!sf_getint("pad1",&pad1)) pad1=1;
    /* padding factor on the first axis */
    if(!sf_getint("gpz",&gpz)) gpz=0;
    /* geophone surface */

    if (adj) { /* migration */
	data = sf_input("in");
	image = sf_output("out");
	sf_settype(image,SF_COMPLEX);

	if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in input");

	if (!sf_histint(data,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(data,"o2",&ox)) ox=0.; 

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* depth samples */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* depth sampling */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putfloat(image,"o1",0.);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n2",nx);
	sf_putfloat(image,"d2",dx);
	sf_putfloat(image,"o2",ox);
	sf_putstring(image,"label2","Distance");
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");
	sf_settype(data,SF_COMPLEX);

	if (!sf_histint(image,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(image,"d1",&dz)) sf_error("No d1= in input");

	if (!sf_histint(image,"n2",&nx))  sf_error("No n2= in input");
	if (!sf_histfloat(image,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(image,"o2",&ox)) ox=0.; 	

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling */

	sf_putint(data,"n1",nt);
	sf_putfloat(data,"d1",dt);
	sf_putfloat(data,"o1",0.);
	sf_putstring(data,"label1","Time");
	sf_putstring(data,"unit1","s");

	sf_putint(data,"n2",nx);
	sf_putfloat(data,"d2",dx);
	sf_putfloat(data,"o2",ox);
	sf_putstring(data,"label2","Distance");
    }

    nz2 = kiss_fft_next_fast_size(nz*pad1);
    nx2 = kiss_fft_next_fast_size(nx);
    nk = nz2*nx2; /*wavenumber*/

    nzx = nz*nx;
    nzx2 = nz2*nx2;
    ntx = nt*nx;
    wfnt = (int)(nt-1)/snap+1;

    if (snap > 0) {
	snaps = sf_output("snaps");
	/* (optional) snapshot file */
	sf_settype(snaps,SF_COMPLEX);
	sf_putint(snaps,"n1",nz);
	sf_putfloat(snaps,"d1",dz);
	sf_putfloat(snaps,"o1",0.);
	sf_putstring(snaps,"label1","Depth");
	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dx);
	sf_putfloat(snaps,"o2",ox);
	sf_putstring(snaps,"label2","Distance");
	sf_putint(snaps,"n3",wfnt);
	sf_putfloat(snaps,"d3",dt*snap);
	sf_putfloat(snaps,"o3",0.);
	sf_putstring(snaps,"label3","Time");
    } else {
	snaps = NULL;
    }

    /* propagator matrices */
    left = sf_input("left");
    right = sf_input("right");
    
    if (!sf_histint(left,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in left",nzx);
    if (!sf_histint(left,"n2",&m2))  sf_error("No n2= in left");
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nk);
    
    lt = sf_complexalloc2(nzx,m2);
    rt = sf_complexalloc2(m2,nk);
    img = sf_complexalloc2(nz,nx);
    dat = sf_complexalloc2(nt,nx);
    if (snap > 0) wvfld = sf_complexalloc3(nz,nx,wfnt);
    else wvfld = NULL;
    geop = (geopar) sf_alloc(1, sizeof(*geop));

    sf_complexread(lt[0],nzx*m2,left);
    sf_complexread(rt[0],m2*nk,right);
    if (adj) sf_complexread(dat[0],ntx,data);
    else sf_complexread(img[0],nzx,image);

    /*close RSF files*/
    sf_fileclose(left);
    sf_fileclose(right);

#ifdef _OPENMP
#pragma omp parallel
{   
    nth = omp_get_num_threads();
}
    sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    if (timer) t0 = gtod_timer();

    /*load constant geopar elements*/
    geop->nx  = nx;
    geop->nz  = nz;
    geop->dx  = dx;
    geop->dz  = dz;
    geop->ox  = ox;
    geop->gpz = gpz;
    geop->nt  = nt;
    geop->dt  = dt;
    geop->snap= snap;
    geop->nzx2= nzx2;
    geop->nk  = nk;
    geop->m2  = m2;
    geop->wfnt= wfnt;

    lrexp(img, dat, adj, lt, rt, geop, pad1, verb, snap, wvfld);

    if (timer)
      {
        t1 = gtod_timer();
        time = t1-t0;
        sf_warning("Time = %lf\n",time);
      }

    if (adj) {
	sf_complexwrite(img[0],nzx,image);
    } else {
	sf_complexwrite(dat[0],ntx,data);
    }
    
    if (snap > 0 && NULL != snaps) {
	sf_complexwrite(wvfld[0][0],wfnt*nx*nz,snaps);
    }

    cfft2_finalize();
    exit(0);
}
