/* 2-D FFT-based zero-offset exploding reflector modeling/migration linear operator */
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

#include "cfft2w.h"
#include "timer.h"
/*******************************************************/
/* wave propagation utils*/
typedef struct Geopar {
    /*geometry parameters*/
    int   nx;
    int   nz;
    float dx;
    float dz;
    int   gpz;
    /*time parameters*/
    int nt;
    float dt;
    int snap;
    int wfnt;
    /*miscellaneous*/
    int pad1;
    bool verb;
    /*propagation matrix*/
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
    /*for tapering */
    int taper;
    float thresh;
} * geopar; /*geometry parameters*/

int lrexp(sf_complex **img, sf_complex **dat, bool adj, sf_complex **lt, sf_complex **rt, geopar geop, sf_complex ***wvfld)
/*< zero-offset exploding reflector modeling/migration >*/
{
    int it, nt, ix, nx, nx2, iz, nz, nz2, nzx2, gpz, wfit=0, snap;
    int im, i, j, m2, ik, nk, pad1;
    float dx, dz;
    sf_complex *curr, **wave, *cwave, *cwavem, c;
    sf_complex *currm;
    bool verb;
    
    /*for tapering*/
    int taper;
    float dkx,dkz,kx0,kz0,kx,kz,ktmp,kx_trs,kz_trs,thresh;
    float *ktp;

    nx  = geop->nx;
    nz  = geop->nz;
    dx  = geop->dx;
    dz  = geop->dz;
    gpz = geop->gpz;
    nt  = geop->nt;
    snap= geop->snap;
    pad1= geop->pad1;
    verb= geop->verb;
    nzx2= geop->nzx2;
    /*nk = geop->nk;*/
    m2  = geop->m2;
    taper = geop->taper;
    thresh = geop->thresh;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    if (nk!=geop->nk) sf_error("nk discrepancy!");

    if (taper) {
        dkz = 1./(nz2*dz); kz0 = -0.5/dz;
        dkx = 1./(nx2*dx); kx0 = -0.5/dx;
        kz_trs = thresh*fabs(kz0);
        kx_trs = thresh*fabs(kx0);
        sf_warning("dkz=%f,dkx=%f,kz0=%f,kx0=%f",dkz,dkx,kz0,kx0);
        sf_warning("nk=%d,nkz=%d,nkx=%d",nk,nz2,nx2);
        sf_warning("Applying kz tapering below %f",kz_trs);
        sf_warning("Applying kx tapering below %f",kx_trs);
        ktp = sf_floatalloc(nk);
        /* constructing the tapering op */
        for (ix=0; ix < nx2; ix++) {
            kx = kx0+ix*dkx;
            for (iz=0; iz < nz2; iz++) {
                kz = kz0+iz*dkz;
                ktmp = 1.;
                if (fabs(kx) > kx_trs)
                    ktmp *= (fabs(kx)>kx_trs)? powf((fabs(kx0)-fabs(kx)+kx_trs)/kx0,2) : 1.;
                if (fabs(kz) > kz_trs)
                    ktmp *= (fabs(kz)>kz_trs)? powf((fabs(kz0)-fabs(kz)+kz_trs)/kz0,2) : 1.;
                ktp[iz+ix*nz2] = ktmp;
            }
        }
    } else {
	ktp = NULL;
    }

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

    if (adj) { /* migration */
        if (snap>0) wfit = (int)(nt-1)/snap; // wfnt-1
	/* time stepping */
	for (it=nt-1; it > -1; it--) {
	    if (verb) sf_warning("it=%d;",it);
	
	    /* 4 - matrix multiplication */
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
		    c = sf_cadd(c,sf_cmul(wave[im][ik],conjf(rt[ik][im]))); //complex multiplies complex
#endif
		}
		cwave[ik] = c;
	    }

	    icfft2(curr,cwave);

            /* 3 - inject data */
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    for (ix=0; ix < nx; ix++) {
#ifdef SF_HAS_COMPLEX_H		
		curr[gpz+ix*nz2] += dat[ix][it];
#else
		curr[gpz+ix*nz2] = sf_cadd(curr[gpz+ix*nz2],dat[ix][it]);
#endif
	    }

            /* 2 - taper */
            if (taper) {
                if (it%taper == 0) {
                    cfft2(curr,cwave);
                    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                        cwave[ik] = cwave[ik]*ktp[ik];
#else
                        cwave[ik] = sf_crmul(cwave[ik],ktp[ik]);
#endif
                    }
                    icfft2(curr,cwave);
                }
            }

            /* 1 - snapshot */
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
		if (snap>0) wfit--;
	    }
	} /*time iteration*/
	/* 0 - generate image */
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		img[ix][iz] = curr[iz+ix*nz2];
	    }
	}
    } else { /* modeling */
	/* 0 - exploding reflector */
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		curr[iz+ix*nz2]=img[ix][iz];
	    }
	}
	if (snap>0) wfit = 0;
	/* time stepping */
	for (it=0; it < nt; it++) {
	    if (verb) sf_warning("it=%d;",it);

            /* 1 - snapshot */
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
		if (snap>0) wfit++;
	    }

            /* 2 - taper */
            if (taper) {
                if (it%taper == 0) {
                    cfft2(curr,cwave);
                    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
                        cwave[ik] = cwave[ik]*ktp[ik];
#else
                        cwave[ik] = sf_crmul(cwave[ik],ktp[ik]);
#endif
                    }
                    icfft2(curr,cwave);
                }
            }

            /* 3 - record data */
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    for (ix=0; ix < nx; ix++) { /* record data */
		dat[ix][it] = curr[gpz+ix*nz2];
	    }

	    /* 4 - matrix multiplication */
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
			c = sf_cadd(c,sf_cmul(lt[im][i], wave[im][j]));
#endif	    
		    }
		    curr[j] = c;
		}
	    }
	}
    }
    if (verb) sf_warning(".");

    cfft2_finalize();
    return 0;
}

int main(int argc, char* argv[])
{
    bool adj,timer,verb;
    int nt, nx, nz, nx2, nz2, nzx, nzx2, ntx, pad1, snap, gpz, wfnt;
    int m2, n2, nk, taper;
#ifdef _OPENMP
    int nth=1;
#endif
    float dt, dx, dz, ox, oz, thresh;
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
    if (!sf_getint("taper",&taper)) taper=0; 
    /* tapering in the frequency domain */
    if (!sf_getfloat("thresh",&thresh)) thresh=0.92; 
    /* tapering threshold */

    if (adj) { /* migration */
	data = sf_input("in");
	image = sf_output("out");
	sf_settype(image,SF_COMPLEX);

	if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in input");
        if (!sf_getfloat("oz",&oz)) oz=0.; 

	if (!sf_histint(data,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(data,"o2",&ox)) ox=0.; 

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* depth samples */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* depth sampling */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putfloat(image,"o1",oz);
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
	if (!sf_histfloat(image,"o1",&oz)) oz=0.;

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

    if (snap > 0) {
        wfnt = (int)(nt-1)/snap+1;
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
        wfnt = 0;
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
    geop->gpz = gpz;
    geop->nt  = nt;
    geop->dt  = dt;
    geop->snap= snap;
    geop->wfnt= wfnt;
    geop->pad1= pad1;
    geop->verb= verb;
    geop->nzx2= nzx2;
    geop->nk  = nk;
    geop->m2  = m2;
    geop->taper= taper;
    geop->thresh= thresh;

    lrexp(img, dat, adj, lt, rt, geop, wvfld);

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

    exit(0);
}
