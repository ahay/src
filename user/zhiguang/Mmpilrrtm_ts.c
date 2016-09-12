/* One-step lowrank RTM with time-shift imaging condition */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/time.h>
#include <math.h>

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

void reflgen(int nzb, int nxb, int spz, int spx,
             int rectz, int rectx, int nrep, /*smoothing parameters*/
			 float *refl/*reflectivity map*/);

static int n1, n2, nkk;
static float wt;

static sf_complex **cc;

#ifdef SF_HAS_FFTW
static fftwf_plan cfg=NULL, icfg=NULL;
#else
static kiss_fft_cfg cfg1, icfg1, cfg2, icfg2;
static kiss_fft_cpx **temp, *ctrace2;
static sf_complex *trace2;
#endif

int cfft2_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny   /* input data size */, 
	       int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{

#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_init_threads();
    sf_warning("Using threaded FFTW3! \n");
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
#endif

#ifndef SF_HAS_FFTW
    int i2;
#endif

    nkk = n1 = kiss_fft_next_fast_size(nx*pad1);
    
#ifndef SF_HAS_FFTW
    cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
    icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);
#endif
  
    n2 = kiss_fft_next_fast_size(ny);

    cc = sf_complexalloc2(n1,n2);
    
#ifndef SF_HAS_FFTW
    cfg2  = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(n2,1,NULL,NULL);
 	
    temp =    (kiss_fft_cpx **) sf_alloc(n2,sizeof(*temp));
    temp[0] = (kiss_fft_cpx *)  sf_alloc(nkk*n2,sizeof(kiss_fft_cpx));
    for (i2=0; i2 < n2; i2++) {
	temp[i2] = temp[0]+i2*nkk;
    }
	
    trace2 = sf_complexalloc(n2);
    ctrace2 = (kiss_fft_cpx *) trace2;
#endif

    *nx2 = n1;
    *ny2 = n2;
	
    wt =  1.0/(n1*n2);
	
    return (nkk*n2);
}

void cfft2(sf_complex *inp /* [n1*n2] */, 
	   sf_complex *out /* [nkk*n2] */)
/*< 2-D FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    if (NULL==cfg) {
      cfg = fftwf_plan_dft_2d(n2,n1,
			      (fftwf_complex *) cc[0], 
			      (fftwf_complex *) out,
			      FFTW_FORWARD, FFTW_MEASURE);
      if (NULL == cfg) sf_error("FFTW failure.");
    }
#endif

    /* FFT centering */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		cc[i2][i1] = ((i2%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1];
#else
		cc[i2][i1] = ((i2%2==0)==(i1%2==0))? inp[i2*n1+i1]:sf_cneg(inp[i2*n1+i1]);
#endif
	}
    }

#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);
#else	
    for (i2=0; i2 < n2; i2++) {
	kiss_fft_stride(cfg1,(kiss_fft_cpx *) cc[i2],temp[i2],1);
    }

    for (i1=0; i1 < nkk; i1++) {
	kiss_fft_stride(cfg2,temp[0]+i1,ctrace2,nkk);
	for (i2=0; i2<n2; i2++) {
	    out[i2*nkk+i1] = trace2[i2];
	}
    }
#endif
}

void icfft2_allocate(sf_complex *inp /* [nkk*n2] */)
/*< allocate inverse transform >*/
{
#ifdef SF_HAS_FFTW
    icfg = fftwf_plan_dft_2d(n2,n1,
			     (fftwf_complex *) inp, 
			     (fftwf_complex *) cc[0],
			     FFTW_BACKWARD, FFTW_MEASURE);
    if (NULL == icfg) sf_error("FFTW failure.");
#endif
}

void icfft2(sf_complex *out /* [n1*n2] */, 
	    sf_complex *inp /* [nkk*n2] */)
/*< 2-D inverse FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    fftwf_execute(icfg);
#else
    for (i1=0; i1 < nkk; i1++) {
	kiss_fft_stride(icfg2,(kiss_fft_cpx *) (inp+i1),ctrace2,nkk);
		
	for (i2=0; i2<n2; i2++) {
	    temp[i2][i1] = ctrace2[i2];
	}
    }
    for (i2=0; i2 < n2; i2++) {
	kiss_fft_stride(icfg1,temp[i2],(kiss_fft_cpx *) cc[i2],1);
    }
#endif
    
    /* FFT centering and normalization*/
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		out[i2*n1+i1] = (((i2%2==0)==(i1%2==0))? wt:-wt) * cc[i2][i1];
#else
		out[i2*n1+i1] = sf_crmul(cc[i2][i1],(((i2%2==0)==(i1%2==0))? wt:-wt));
#endif
	}
    }
}

void cfft2_finalize()
/*< clean up fftw >*/
{
/* make sure everything is back to its pristine state */
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_cleanup_threads();
#endif
    fftwf_destroy_plan(cfg);
    fftwf_destroy_plan(icfg);
    fftwf_cleanup();
    cfg=NULL;
    icfg=NULL;
#else
    free(cfg1); cfg1=NULL;
    free(icfg1); icfg1=NULL;
    free(cfg2); cfg2=NULL;
    free(icfg2); icfg2=NULL;
#endif

    free(*cc);
    free(cc);
}

int main(int argc, char *argv[])
{
	bool wantwf, verb;

	int ix, iz, is, it, wfit, im, ik, i, j, itau;
    int ns, nx, nz, nt, wfnt, rnx, rnz, nzx, rnzx, vnx, ntau, htau, nds;
	int scalet, snap, snapshot, fnx, fnz, fnzx, nk, nb;
	int rectx, rectz, gpz, n, m, pad1, trunc, spx, spz;

	float dt, t0, z0, dz, x0, dx, s0, ds, wfdt, srctrunc;
    float dtau, tau0, tau;

	char *path1, *path2, number[5], *left, *right;

	double tstart, tend;
	struct timeval tim;
	
	sf_complex c, **lt, **rt;
	sf_complex *ww, **dd;
	float ***img1, **img2, ***mig1, **mig2;
    float *rr, **ccr, **sill, ***fwf, ***bwf;
	sf_complex *cwave, *cwavem, **wave, *curr;

	sf_axis at, ax, az, atau;

	sf_file Fdat, Fsrc, Fimg1, Fimg2;
	sf_file Ffwf, Fbwf, Fvel;
	sf_file Fleft, Fright;

	int cpuid, numprocs, nth;
    float *sendbuf, *recvbuf;
	MPI_Comm comm=MPI_COMM_WORLD;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &cpuid);
	MPI_Comm_size(comm, &numprocs);

	sf_init(argc, argv);

#ifdef _OPENMP
#pragma omp parallel
	{
		nth=omp_get_num_threads();
	}
	sf_warning(">>> Using %d threads <<<", nth);
#endif

	gettimeofday(&tim, NULL);
	tstart=tim.tv_sec+(tim.tv_usec/1000000.0);

	if(!sf_getbool("wantwf", &wantwf)) wantwf=false;
    if(!sf_getbool("verb", &verb)) verb=false;
	if(!sf_getint("pad1", &pad1)) pad1=1;
	/* padding factor on the first axis */

	if(!sf_getint("nb", &nb)) sf_error("Need nb= ");
	if(!sf_getfloat("srctrunc", &srctrunc)) srctrunc=0.4;
	if(!sf_getint("rectx", &rectx)) rectx=2;
	if(!sf_getint("rectz", &rectz)) rectz=2;

	if(!sf_getint("scalet", &scalet)) scalet=1;
	if(!sf_getint("snap", &snap)) snap=100;
	/* interval of the output wavefield */
	if(!sf_getint("snapshot", &snapshot)) snapshot=0;
	/* print out the wavefield snapshots of this shot */
    if(!sf_getint("nds", &nds)) sf_error("Need nds=!");
    
    /* source and receiver positions */
	if(!sf_getint("gpz", &gpz)) sf_error("Need gpz=");
	if(!sf_getint("spx", &spx)) sf_error("Need spx=");
	if(!sf_getint("spz", &spz)) sf_error("Need spz=");
    
    /* tau parameters */
    if(!sf_getint("ntau", &ntau)) sf_error("Need ntau=");
    if(!sf_getfloat("dtau", &dtau)) sf_error("Need dtau=");
    if(!sf_getfloat("tau0", &tau0)) sf_error("Need tau0=");

	/* input/output files */
	Fdat=sf_input("--input");
	Fimg1=sf_output("--output");
    Fimg2=sf_output("Fimg2");
    Fsrc=sf_input("Fsrc");
    Fvel=sf_input("Fpadvel");

	if(wantwf){
		Ffwf=sf_output("Ffwf");
        Fbwf=sf_output("Fbwf");
	}

	at=sf_iaxa(Fsrc, 1); nt=sf_n(at); dt=sf_d(at); t0=sf_o(at);
    ax=sf_iaxa(Fvel, 2); vnx=sf_n(ax); x0=sf_o(ax);
	az=sf_iaxa(Fvel, 1); rnz=sf_n(az); dz=sf_d(az); z0=sf_o(az);
    if(!sf_histint(Fdat, "n2", &rnx)) sf_error("Need n2= in input!");
    if(!sf_histfloat(Fdat, "d2", &dx)) sf_error("Need d2= in input!");
    if(!sf_histint(Fdat, "n3", &ns)) sf_error("Need n3= in input!");
    if(!sf_histfloat(Fdat, "d3", &ds)) sf_error("Need d3= in input!");
    if(!sf_histfloat(Fdat, "o3", &s0)) sf_error("Need o3= in input!");
    
    wfnt=(nt-1)/scalet+1;
    wfdt=dt*scalet;
    
    /* double check the geometry parameters */
    if(nds != (int)(ds/dx)) sf_error("Need ds/dx= %d", nds);
	sf_warning("s0=%g, x0+(rnx-1)*dx/2=%g", s0, x0+(rnx-1)*dx/2);
    //if(s0 != x0+(rnx-1)*dx/2) sf_error("Wrong origin information!");
    if(vnx != nds*(ns-1)+rnx) sf_error("Wrong dimension in x axis!");

    /* set up the output files */
    atau=sf_iaxa(Fsrc, 1);
    sf_setn(atau, ntau);
    sf_setd(atau, dtau);
    sf_seto(atau, tau0);
    sf_setlabel(atau, "Tau");
    sf_setunit(atau, "s");
    
    sf_oaxa(Fimg1, az, 1);
    sf_oaxa(Fimg1, ax, 2);
    sf_oaxa(Fimg1, atau, 3);
    sf_oaxa(Fimg2, az, 1);
    sf_oaxa(Fimg2, ax, 2);
    sf_putint(Fimg2, "n3", 1);
    sf_settype(Fimg1, SF_FLOAT);
    sf_settype(Fimg2, SF_FLOAT);
    
    if(wantwf){
		sf_setn(ax, rnx);
        sf_seto(ax, -(rnx-1)*dx/2.0);
        sf_oaxa(Ffwf, az, 1);
        sf_oaxa(Ffwf, ax, 2);
        sf_putint(Ffwf, "n3", (wfnt-1)/snap+1);
        sf_putfloat(Ffwf, "d3", snap*wfdt);
        sf_putfloat(Ffwf, "o3", t0);
        sf_putstring(Ffwf, "label3", "Time");
        sf_putstring(Ffwf, "unit3", "s");
        sf_settype(Ffwf, SF_FLOAT);
        
        sf_oaxa(Fbwf, az, 1);
        sf_oaxa(Fbwf, ax, 2);
        sf_putint(Fbwf, "n3", (wfnt-1)/snap+1);
        sf_putfloat(Fbwf, "d3", -snap*wfdt);
        sf_putfloat(Fbwf, "o3", (wfnt-1)*wfdt);
        sf_putstring(Fbwf, "label3", "Time");
        sf_putstring(Fbwf, "unit3", "s");
        sf_settype(Fbwf, SF_FLOAT);
	}
	
    nx=rnx+2*nb; nz=rnz+2*nb;
	nzx=nx*nz; rnzx=rnz*rnx;
    nk=cfft2_init(pad1, nz, nx, &fnz, &fnx);
	fnzx=fnz*fnx;
    
	/* print axies parameters for double check */
    sf_warning("cpuid=%d, numprocs=%d", cpuid, numprocs);
	sf_warning("nt=%d, dt=%g, scalet=%d, wfnt=%d, wfdt=%g",nt, dt, scalet, wfnt, wfdt);
	sf_warning("vnx=%d, nx=%d, dx=%g, nb=%d, rnx=%d", vnx, nx, dx, nb, rnx);
	sf_warning("nz=%d, rnz=%d, dz=%g, z0=%g", nz, rnz, dz, z0);
	sf_warning("spx=%d, spz=%d, gpz=%d", spx, spz, gpz);
	sf_warning("ns=%d, ds=%g, s0=%g", ns, ds, s0);
    sf_warning("ntau=%d, dtau=%g, tau0=%g", ntau, dtau, tau0);
    sf_warning("nzx=%d, fnzx=%d, nk=%d", nzx, fnzx, nk);

	/* allocate storage and read data */
	ww=sf_complexalloc(nt);
	sf_complexread(ww, nt, Fsrc);
	sf_fileclose(Fsrc);
	
    gpz=gpz+nb;
    spz=spz+nb;
    spx=spx+nb;
    trunc=srctrunc/dt+0.5;
    
	dd=sf_complexalloc2(nt, rnx);
	rr=sf_floatalloc(nzx);
	reflgen(nz, nx, spz, spx, rectz, rectx, 0, rr);
    
    fwf=sf_floatalloc3(rnz, rnx, wfnt);
    bwf=sf_floatalloc3(rnz, rnx, wfnt);
    img1=sf_floatalloc3(rnz, vnx, ntau);
    img2=sf_floatalloc2(rnz, vnx);
    mig1=sf_floatalloc3(rnz, rnx, ntau);
    mig2=sf_floatalloc2(rnz, rnx);
    
    ccr=sf_floatalloc2(rnz, rnx);
    sill=sf_floatalloc2(rnz, rnx);
    
    curr=sf_complexalloc(fnzx);
	cwave=sf_complexalloc(nk);
	cwavem=sf_complexalloc(nk);
    icfft2_allocate(cwavem);
    
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, itau)
#endif
    for(ix=0; ix<vnx; ix++){
        for(iz=0; iz<rnz; iz++){
            img2[ix][iz]=0.;
            for(itau=0; itau<ntau; itau++){
                img1[itau][ix][iz]=0.;
            }
        }
    }

	path1=sf_getstring("path1");
	path2=sf_getstring("path2");
	if(path1==NULL) path1="./mat/left";
	if(path2==NULL) path2="./mat/right";
	/* shot loop */
	for (is=cpuid; is<ns; is+=numprocs){
		/* construct the names of left and right matrices */
		left=sf_charalloc(strlen(path1));
		right=sf_charalloc(strlen(path2));
		strcpy(left, path1);
		strcpy(right, path2);
		sprintf(number, "%d", is+1);
		strcat(left, number);
		strcat(right, number);

		Fleft=sf_input(left);
		Fright=sf_input(right);
		
		if(!sf_histint(Fleft, "n1", &n) || n != nzx) sf_error("Need n1=%d in Fleft", nzx);
		if(!sf_histint(Fleft, "n2", &m)) sf_error("No n2 in Fleft");
		if(!sf_histint(Fright, "n1", &n) || n != m) sf_error("Need n1=%d in Fright", m);
		if(!sf_histint(Fright, "n2", &n) || n != nk) sf_error("Need n2=%d in Fright", nk);
		
		/* allocate storage for each shot migration */
		lt=sf_complexalloc2(nzx, m);
		rt=sf_complexalloc2(m, nk);
		sf_complexread(lt[0], nzx*m, Fleft);
		sf_complexread(rt[0], m*nk, Fright);
        sf_fileclose(Fleft);
		sf_fileclose(Fright);
        
        /* read data */
        sf_seek(Fdat, is*rnx*nt*sizeof(float complex), SEEK_SET);
        sf_complexread(dd[0], rnx*nt, Fdat);
        
        /* initialize curr and imaging variables */
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
		for(iz=0; iz<fnzx; iz++){
			curr[iz]=sf_cmplx(0.,0.);
		}
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, itau)
#endif
        for(ix=0; ix<rnx; ix++){
            for(iz=0; iz<rnz; iz++){
                mig2[ix][iz]=0.;
                ccr[ix][iz]=0.;
                sill[ix][iz]=0.;
                for(itau=0; itau<ntau; itau++){
                    mig1[itau][ix][iz]=0.;
                }
            }
        }
        
        /* wave */
		wave=sf_complexalloc2(fnzx, m);
        
        /* snapshot */
        if(wantwf && is== snapshot) wantwf=true;
        else wantwf=false;
		
		wfit=0;
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation it=%d/%d",it+1, nt);
			
			cfft2(curr, cwave);
			for(im=0; im<m; im++){
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
				for(ik=0; ik<nk; ik++){
#ifdef SF_HAS_COMPLEX_H
					cwavem[ik]=cwave[ik]*rt[ik][im];
#else
					cwavem[ik]=sf_cmul(cwave[ik],rt[ik][im]);
#endif
				}
				icfft2(wave[im],cwavem);
			}

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, im, c) shared(curr, it)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					i=iz+ix*nz;
					j=iz+ix*fnz;
					
					if(it<trunc){
#ifdef SF_HAS_COMPLEX_H
						c=ww[it]*rr[i];
#else
						c=sf_crmul(ww[it],rr[i]);
#endif
					}else{
						c=sf_cmplx(0.,0.);
					}
					
//                    c += curr[j];
                    
					for(im=0; im<m; im++){
#ifdef SF_HAS_COMPLEX_H
						c += lt[im][i]*wave[im][j];
#else
						c += sf_cmul(lt[im][i], wave[im][j]);
#endif
					}
					curr[j]=c;
				}
			}
			
			if(it%scalet==0){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
				for(ix=0; ix<rnx; ix++){
                    for(iz=0; iz<rnz; iz++){
                        fwf[wfit][ix][iz]=crealf(curr[(ix+nb)*fnz+(iz+nb)]);
                    }
                }
                wfit++;
            }
        } //end of it
        
        /* check wfnt */
        if(wfit != wfnt) sf_error("At this point, wfit should be equal to wfnt");

        
        /* backward propagation starts from here... */
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
		for(iz=0; iz<fnzx; iz++){
			curr[iz]=sf_cmplx(0.,0.);
		}
        
        wfit=wfnt-1;
        for(it=nt-1; it>=0; it--){
            if(verb) sf_warning("Backward propagation it=%d/%d",it+1, nt);
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
            for(ix=0; ix<rnx; ix++){
                curr[(ix+nb)*fnz+gpz]+=dd[ix][it];
            }
            
            cfft2(curr, cwave);
            
			for(im=0; im<m; im++){
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
				for(ik=0; ik<nk; ik++){
#ifdef SF_HAS_COMPLEX_H
					cwavem[ik]=cwave[ik]*conjf(rt[ik][im]);
#else
					cwavem[ik]=sf_cmul(cwave[ik],conjf(rt[ik][im]));
#endif
				}
				icfft2(wave[im],cwavem);
			}
            
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, im, c) shared(curr, it)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					i=iz+ix*nz;
					j=iz+ix*fnz;
					
//                    c=curr[j];
                      c=sf_cmplx(0.,0.);
					
					for(im=0; im<m; im++){
#ifdef SF_HAS_COMPLEX_H
						c += conjf(lt[im][i])*wave[im][j];
#else
						c += sf_cmul(conjf(lt[im][i]), wave[im][j]);
#endif
					}
					curr[j]=c;
				}
			}
			
			if(it%scalet==0){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
				for(ix=0; ix<rnx; ix++){
                    for(iz=0; iz<rnz; iz++){
                        bwf[wfit][ix][iz]=crealf(curr[(ix+nb)*fnz+(iz+nb)]);
                        ccr[ix][iz] += fwf[wfit][ix][iz]*bwf[wfit][ix][iz];
                        sill[ix][iz] += fwf[wfit][ix][iz]*fwf[wfit][ix][iz];
                    }
                }
                wfit--;
            }
        } //end of it
        if(wfit != -1) sf_error("Check program! The final wfit should be -1!");

        /* free storage */
        free(*rt); free(rt);
        free(*lt); free(lt);
        free(*wave); free(wave);
        free(left); free(right);
        
        /* normalized image */
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
        for (ix=0; ix<rnx; ix++){
            for(iz=0; iz<rnz; iz++){
                mig2[ix][iz]=ccr[ix][iz]/(sill[ix][iz]+SF_EPS);
		//		sill[ix][iz]=0.;
            }
        }
        
        /* calculate time shift gathers */
        for(itau=0; itau<ntau; itau++){
			sf_warning("itau/ntau=%d/%d", itau+1, ntau);
            tau=itau*dtau+tau0;
            htau=tau/wfdt;

            for(it=abs(htau); it<wfnt-abs(htau); it++){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
                for(ix=0; ix<rnx; ix++){
                    for(iz=0; iz<rnz; iz++){
                        mig1[itau][ix][iz]+=fwf[it+htau][ix][iz]*bwf[it-htau][ix][iz];
		//				sill[ix][iz]+=fwf[it+htau][ix][iz]*fwf[it+htau][ix][iz];
                    } // end of iz
                } // end of ix
            } // end of it
			
/*
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif */
			/* source illumination
			for(ix=0; ix<rnx; ix++){
				for(iz=0; iz<rnz; iz++){
					mig1[itau][ix][iz] = mig1[itau][ix][iz]/(sill[ix][iz]+SF_EPS);
				}
			} */
        } //end of itau
        
        /* output wavefield snapshot */
        if(wantwf){
            for(it=0; it<wfnt; it++){
                if(it%snap==0){
                    sf_floatwrite(fwf[it][0], rnzx, Ffwf);
                    sf_floatwrite(bwf[wfnt-1-it][0], rnzx, Fbwf);
                }
            }
            sf_fileclose(Ffwf);
            sf_fileclose(Fbwf);
        }
        
        /* add all the shot images that are on the same node */
#ifdef _OPENMP
#pragma omp parallel for private(itau, ix, iz)
#endif
        for(itau=0; itau<ntau; itau++){
            for(ix=0; ix<rnx; ix++){
                for(iz=0; iz<rnz; iz++){
                    img1[itau][ix+is*nds][iz] += mig1[itau][ix][iz];
                }
            }
        }
        
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
        for(ix=0; ix<rnx; ix++){
            for(iz=0; iz<rnz; iz++){
                img2[ix+is*nds][iz] += mig2[ix][iz];
            }
        }
	}
    ////////////////end of ishot
	MPI_Barrier(comm);
    
    cfft2_finalize();
    sf_fileclose(Fdat);
    
    free(ww); free(rr);
	free(*dd); free(dd);
	free(cwave); free(cwavem); free(curr);
    free(*ccr); free(ccr);
    free(*sill); free(sill);
    free(**fwf); free(*fwf); free(fwf);
    free(**bwf); free(*bwf); free(bwf);
    free(**mig1); free(*mig1); free(mig1);
    free(*mig2); free(mig2);
    
    /* sum image */
    if(cpuid==0){
        sendbuf=(float *)MPI_IN_PLACE;
        recvbuf=img1[0][0];
    }else{
        sendbuf=img1[0][0];
        recvbuf=NULL;
    }
    MPI_Reduce(sendbuf, recvbuf, ntau*vnx*rnz, MPI_FLOAT, MPI_SUM, 0, comm);
    
    if(cpuid==0){
        sendbuf=MPI_IN_PLACE;
        recvbuf=img2[0];
    }else{
        sendbuf=img2[0];
        recvbuf=NULL;
    }
    MPI_Reduce(sendbuf, recvbuf, vnx*rnz, MPI_FLOAT, MPI_SUM, 0, comm);
    
    /* output image */
    if(cpuid==0){
        sf_floatwrite(img1[0][0], ntau*vnx*rnz, Fimg1);
        sf_floatwrite(img2[0], vnx*rnz, Fimg2);
    }
	MPI_Barrier(comm);

	sf_fileclose(Fimg1);
    sf_fileclose(Fimg2);
    free(**img1); free(*img1); free(img1);
    free(*img2); free(img2);
    
	gettimeofday(&tim, NULL);
	tend=tim.tv_sec+(tim.tv_usec/1000000.0);
	sf_warning(">> The computing time is %.3lf minutes <<", (tend-tstart)/60.);

	MPI_Finalize();
	exit(0);
}

void reflgen(int nzb, int nxb, int spz, int spx,
             int rectz, int rectx, int nrep, /*smoothing parameters*/
	     float *refl/*reflectivity map*/)
/*< Generate reflectivity map with smoothing >*/
{   
    int i, j, i0, irep;
    int nzx=nzb*nxb;
    sf_triangle tr;
    int n[2],s[2],rect[2];
    bool diff[2],box[2];

    n[0]=nzb; n[1]=nxb;
    s[0]=1;   s[1]=nzb;
    rect[0]=rectz; rect[1]=rectx;
    diff[0]=false; diff[1]=false;
    box[0]=false; box[1]=false;
    
    j=spx*nzb+spz; /*point source position*/
    refl[j]=1;
    
    /* 2-d triangle smoothing */
    for (i=0;i<2;i++) {
      if (rect[i] <= 1) continue;
      tr = sf_triangle_init (rect[i],n[i],box[i]);
      for (j=0; j < nzx/n[i]; j++) {
	i0 = sf_first_index (i,j,2,n,s);
	for (irep=0; irep < nrep; irep++) {
	  sf_smooth2 (tr,i0,s[i],diff[i],refl);
	}
      }
      sf_triangle_close(tr);
    }
}
