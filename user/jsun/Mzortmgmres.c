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

#include "cfft2nsps.h"
#include "timer.h"
#include "zolsrtm2.h"
#include "cgmres.h"

int main(int argc, char* argv[])
{
    bool adj,timer,verb,gmres;
    int nt, nx, nz, nx2, nz2, nzx, nzx2, ntx, pad1, snap, gpz, wfnt, i;
    int m2, n2, nk, nth=1;
    int niter, mem;
    float dt, dx, dz, ox;
    sf_complex *img, *imgout, *dat, **lt1, **rt1, **lt2, **rt2, ***wvfld;
    sf_file data, image, leftf, rightf, leftb, rightb, snaps;
    double time=0.,t0=0.,t1=0.;
    geopar geop;

    sf_init(argc,argv);

    /* essentially doing imaging */
    adj = true;

    if (!sf_getbool("gmres",&gmres)) gmres=false;
    if (gmres) {
      if (!sf_getint("niter",&niter)) niter=10;
      if (!sf_getint("mem",&mem)) mem=20;
    }

    if(! sf_getbool("timer",&timer)) timer=false;
    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */
    if (!sf_getint("pad1",&pad1)) pad1=1;
    /* padding factor on the first axis */
    if(!sf_getint("gpz",&gpz)) gpz=0;
    /* geophone surface */

    /* adj */
    if (!sf_getint("nz",&nz)) sf_error("Need nz=");
    /* depth samples */
    if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
    /* depth sampling */

    /* for */
    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* time sampling */

    if (adj) { /* migration */
	data = sf_input("in");
	image = sf_output("out");
	sf_settype(image,SF_COMPLEX);

	if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in input");

	if (!sf_histint(data,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(data,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(data,"o2",&ox)) ox=0.; 

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
    leftf = sf_input("leftf");
    rightf = sf_input("rightf");
    
    if (!sf_histint(leftf,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in leftf",nzx);
    if (!sf_histint(leftf,"n2",&m2))  sf_error("No n2= in leftf");
    if (!sf_histint(rightf,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in rightf",m2);
    if (!sf_histint(rightf,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in rightf",nk);

    leftb = sf_input("leftb");
    rightb = sf_input("rightb");
    
    if (!sf_histint(leftb,"n1",&n2) || n2 != nzx) sf_error("Need n1=%d in leftb",nzx);
    if (!sf_histint(leftb,"n2",&m2))  sf_error("No n2= in leftb");
    if (!sf_histint(rightb,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in rightb",m2);
    if (!sf_histint(rightb,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in rightb",nk);
    
    lt1 = sf_complexalloc2(nzx,m2); /* propagator for forward modeling */
    rt1 = sf_complexalloc2(m2,nk);
    lt2 = sf_complexalloc2(nzx,m2); /* propagator for backward imaging */
    rt2 = sf_complexalloc2(m2,nk);
    img = sf_complexalloc(nz*nx);
    dat = sf_complexalloc(nt*nx);
    imgout = sf_complexalloc(nz*nx);
    if (snap > 0) wvfld = sf_complexalloc3(nz,nx,wfnt);
    else wvfld = NULL;
    geop = (geopar) sf_alloc(1, sizeof(*geop));

    sf_complexread(lt1[0],nzx*m2,leftf);
    sf_complexread(rt1[0],m2*nk,rightf);
    sf_complexread(lt2[0],nzx*m2,leftb);
    sf_complexread(rt2[0],m2*nk,rightb);
    if (adj) sf_complexread(dat,ntx,data);
    else sf_complexread(img,nzx,image);

    /*close RSF files*/
    sf_fileclose(leftf);
    sf_fileclose(rightf);
    sf_fileclose(leftb);
    sf_fileclose(rightb);

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
    geop->pad1= pad1;
    geop->verb= verb;

    /* first get the Q-compensated image: B[d] */
    lrexp(img, dat, adj, lt2, rt2, geop, wvfld);

    /* performing the least-squares optimization: ||{BF}[m] - B[d]|| */
    if (gmres) {
      /* disabling snapshots */
      geop->snap= 0;
      lrexp_init(lt1,rt1,lt2,rt2);

      sf_warning(">>>>>> Using GMRES(m) <<<<<<");
      cgmres_init(nzx,mem);
      cgmres( img, imgout, lrexp_op, geop, niter, 0.01*SF_EPS, true);
      /*lrexp_op(nzx, img, imgout, geop);*/
      lrexp_close();
    } else {
      for (i=0; i<nzx; i++)
	imgout[i] = img[i];
    }

    if (timer)
      {
        t1 = gtod_timer();
        time = t1-t0;
        sf_warning("Time = %lf\n",time);
      }

    if (adj) {
	sf_complexwrite(imgout,nzx,image);
    } else {
	sf_complexwrite(dat,ntx,data);
    }
    
    if (snap > 0 && NULL != snaps) {
	sf_complexwrite(wvfld[0][0],wfnt*nx*nz,snaps);
    }

    exit(0);
}
