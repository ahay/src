/* Conventional FWI misfit and gradient calculation using one-step low-rank wave extrapolation */
/*
 Copyright (C) 2016 University of Texas at Austin
 
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

#include "fftomp.h"

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

int main(int argc, char* argv[])
{
	bool verb, adjtest, scomp;
	int function;

	int ix, iz, is, ir, it, wfit, im, ik, i, j;
	int ns, nr, nt, wfnt, nx, nz, nzx, rnx, rnz, rnzx;
	int scalet, fnx, fnz, fny, fnzx, nk, nb;
	int rectx, rectz, n, m, m2, pad1, ds_v, s0_v, sx, sz, dr_v, r0_v, rz;
	int cpuid, numprocs;

	float dt, t0, dr, r0, ds, s0, z0, dz, x0, dx, wfdt;
	float w0, fcost=0.,gamma;
	float kx,kz,k2,*kk1, *kk2;
	float kz0, kx0, dkz, dkx;
	float *rr, *vv, *qq, *grad;
	float *sendbuf, *recvbuf;
	float a1, a2, v, q;

	sf_complex wf, forw, adjt;
	sf_complex c, **lt, **rt, **ltb, **rtb;
	sf_complex *ww, **pp, **dd;
	sf_complex *cwave, *cwavem, **wave, **wave2, *curr, ***wfl;
	sf_complex *cfrac, *cfrac1, *cfrac2;

	sf_axis at, ax, az;

	sf_file Fsrc, Fvel, Fq, Fleft, Fright, Fleftb, Frightb, Fdat, Fmisfit, Fgrad;
	sf_file Fwav, Fwav2, Fres;
	FILE *swap;
	
	MPI_Comm comm=MPI_COMM_WORLD;

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	sf_init(argc, argv);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &cpuid);
	MPI_Comm_size(comm, &numprocs);

	if(!sf_getint("function", &function)) function=3;
	/* if 1, forward modeling; if 2, only calculate misfit; if 3, calculate gradient */

	Fvel=sf_input("Fvel"); /* velocity model */
	Fsrc=sf_input("Fsrc"); /* wavelet */
	Fleft=sf_input("Fleft"); /* left matrix */
	Fright=sf_input("Fright"); /* left matrix */
	if(function==1){
		swap=fopen("tempswap.bin","wb+");
		Fdat=sf_output("Fdat");
	}else{
		Fdat=sf_input("Fdat");
		Fmisfit=sf_output("Fmisfit");
	}
	if(function==3){
		Fleftb=sf_input("Fleftb");
		Frightb=sf_input("Frightb");
		Fq=sf_input("Fq");
		Fgrad=sf_output("Fgrad");
	}
  
	Fwav=sf_output("Fwav");
	Fwav2=sf_output("Fwav2");
	Fres=sf_output("Fres");
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	ax=sf_iaxa(Fvel, 2); rnx=sf_n(ax); dx=sf_d(ax); x0=sf_o(ax);
	az=sf_iaxa(Fvel, 1); rnz=sf_n(az); dz=sf_d(az); z0=sf_o(az);
	at=sf_iaxa(Fsrc, 1); nt=sf_n(at); dt=sf_d(at); t0=sf_o(at);

	if(!sf_getint("ns", &ns)) sf_error("shot number required"); /* shot number */
	if(!sf_getfloat("ds", &ds)) sf_error("shot interval required"); /* shot interval */
	if(!sf_getfloat("s0", &s0)) sf_error("shot origin required"); /* shot origin */
	if(!sf_getint("sz", &sz)) sz=5; /* source depth */
	if(!sf_getint("nr", &nr)) nr=rnx; /* number of receiver */
	if(!sf_getfloat("dr", &dr)) dr=dx; /* receiver interval */
	if(!sf_getfloat("r0", &r0)) r0=x0; /* receiver origin */
	if(!sf_getint("rz", &rz)) rz=5; /* receiver depth */

	if(!sf_getbool("verb", &verb)) verb=false; /* verbosity flag */
	if(!sf_getbool("scomp", &scomp)) scomp=false; /* source wavefield compensation flag */
	if(!sf_getbool("adjtest", &adjtest)) adjtest=false; /* test of adjointness */
	if(!sf_getint("pad1", &pad1)) pad1=1;
	/* padding factor on the first axis */
	if(!sf_getint("nb", &nb)) sf_error("Need nb=");
	if(!sf_getfloat("w0", &w0)) sf_error("Need w0=");
	/* reference frequency */
	if(!sf_getint("rectx", &rectx)) rectx=3;
	if(!sf_getint("rectz", &rectz)) rectz=3;
	if(!sf_getint("scalet", &scalet)) scalet=1;
	/* time interval */

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	nx=rnx+2*nb; nz=rnz+2*nb;
	nzx=nx*nz; rnzx=rnz*rnx;
	
	wfnt=(nt-1)/scalet+1;
	wfdt=dt*scalet;

	ds_v=ds/dx+0.5;
	s0_v=s0/dx+0.5+nb;
	sz += nb;
	
	dr_v=dr/dx+0.5;
	r0_v=r0/dx+0.5+nb;
	rz += nb;

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	if(function==1){
		sf_setlabel(at, "Time");
		sf_setunit(at, "s");
		sf_oaxa(Fdat, at, 1);

		sf_putint(Fdat, "n2", nr);
		sf_putfloat(Fdat, "d2", dr);
		sf_putfloat(Fdat, "o2", r0);
		sf_putstring(Fdat, "label2", "Receiver");
		sf_putstring(Fdat, "unit2", "km");

		sf_putint(Fdat, "n3", ns);
		sf_putfloat(Fdat, "d3", ds);
		sf_putfloat(Fdat, "o3", s0);
		sf_putstring(Fdat, "label3", "Shot");
		sf_putstring(Fdat, "unit3", "km");
		
		sf_settype(Fdat, SF_COMPLEX);
	}else{
		sf_putint(Fmisfit, "n1", 1);
		sf_putint(Fmisfit, "n2", 1);
		sf_settype(Fmisfit, SF_FLOAT);
	}

	if(function==3){
		sf_oaxa(Fgrad, az, 1);
		sf_oaxa(Fgrad, ax, 2);
		sf_putint(Fgrad, "n3", 1);
		sf_settype(Fgrad, SF_FLOAT);
	}

	sf_putint(Fwav, "n1", rnz);
	sf_putint(Fwav, "n2", rnx);
	sf_putint(Fwav, "n3", (nt-1)/50+1);
	sf_settype(Fwav, SF_COMPLEX);

	sf_putint(Fwav2, "n1", rnz);
	sf_putint(Fwav2, "n2", rnx);
	sf_putint(Fwav2, "n3", (nt-1)/50+1);
	sf_settype(Fwav2, SF_COMPLEX);

	sf_putint(Fres, "n1", nt);
	sf_putint(Fres, "n2", nr);
	sf_putint(Fres, "n3", 1);
	sf_settype(Fres, SF_COMPLEX);

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	nk=fft_init(pad1, nz, nx, 1, &fnz, &fnx, &fny, false, true);
	fnzx=fnz*fnx;

	dkz = 1./(fnz*dz); kz0 = -0.5/dz;
	dkx = 1./(fnx*dx); kx0 = -0.5/dx;

	if(!sf_histint(Fleft, "n1", &n) || n!=nzx) sf_error("Need n1=%d in Fleft", nzx);
	if(!sf_histint(Fleft, "n2", &m)) sf_error("No n2= in Fleft");
	if(!sf_histint(Fright, "n1", &n) || n!=m) sf_error("Need n1=%d in Fright", m);
	if(!sf_histint(Fright, "n2", &n) || n!=nk) sf_error("Need n2=%d in Fright", nk);

	lt=sf_complexalloc2(nzx, m);
	rt=sf_complexalloc2(m, nk);
	sf_complexread(lt[0], nzx*m, Fleft);
	sf_complexread(rt[0], m*nk, Fright);
	sf_fileclose(Fleft);
	sf_fileclose(Fright);

	ww=sf_complexalloc(nt);
	sf_complexread(ww, nt, Fsrc);
	sf_fileclose(Fsrc);

	rr=sf_floatalloc(nzx);
	pp=sf_complexalloc2(nt, nr);
	if(function != 1) dd=sf_complexalloc2(nt, nr);
	if(function==3){
		vv=sf_floatalloc(rnzx);
		qq=sf_floatalloc(rnzx);
		sf_floatread(vv, rnzx, Fvel);
		sf_floatread(qq, rnzx, Fq);
		sf_fileclose(Fvel);
		sf_fileclose(Fq);

		if(!sf_histint(Fleftb, "n2", &m2)) sf_error("No n2= in Fleftb");
		ltb=sf_complexalloc2(nzx, m2);
		rtb=sf_complexalloc2(m2, nk);
		sf_complexread(ltb[0], nzx*m2, Fleftb);
		sf_complexread(rtb[0], m2*nk, Frightb);
		sf_fileclose(Fleftb);
		sf_fileclose(Frightb);

		kk1=sf_floatalloc(fnzx);
		kk2=sf_floatalloc(fnzx);
		cfrac =sf_complexalloc(nk);
		cfrac1=sf_complexalloc(nk);
		cfrac2=sf_complexalloc(nk);

		wfl=sf_complexalloc3(rnz, rnx, wfnt);
		grad=sf_floatalloc(rnzx);
		memset(grad, 0., rnzx*sizeof(float));

		gamma = 0.;
		for(ix=0; ix<rnzx; ix++){
			qq[ix]=atanf(1./qq[ix])/SF_PI;
			gamma+=qq[ix];
		}
		gamma/=rnzx;

        /* calculate fractional Lapl */ 
		for (ix=0; ix<fnx; ix++) {
			kx=2.*SF_PI*(kx0+ix*dkx);
			for (iz=0; iz<fnz; iz++) {
				kz=2.*SF_PI*(kz0+iz*dkz);
				j=ix*fnz+iz;
				k2=kx*kx+kz*kz;
				kk1[j]=powf(k2,gamma+0.5);
				kk2[j]=powf(k2,0.5*gamma+0.5);
			}
		}
		
		wave2=sf_complexalloc2(fnzx, m2);
	} // function==3

	curr=sf_complexalloc(fnzx);
	cwave=sf_complexalloc(nk);
	cwavem=sf_complexalloc(nk);
	wave=sf_complexalloc2(fnzx, m);

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/* shot loop */
	for(is=cpuid; is<ns; is+=numprocs){
		if(function==1) sf_warning("########### is=%d, ns=%d ############", is+1, ns);

		sx=s0_v+is*ds_v;
		reflgen(nz, nx, sz, sx, rectz, rectx, 0., rr);

		if(adjtest && is==ns/2) {
			forw=sf_cmplx(0.,0.);
			adjt=sf_cmplx(0.,0.);
		}

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		/* data misfit */
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
		for(iz=0; iz<fnzx; iz++)
			curr[iz]=sf_cmplx(0.,0.);
		
		wfit=0;
		for(it=0; it<nt; it++){
			if(verb) sf_warning("Forward propagation it=%d/%d", it+1, nt);

			/* forward Fourier transform */
			fft(curr, cwave);

			/* right matrix */
			for(im=0; im<m; im++){
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif 
				for(ik=0; ik<nk; ik++){
#ifdef SF_HAS_COMPLEX_H
					cwavem[ik]=cwave[ik]*rt[ik][im];
#else
					cwavem[ik]=sf_cmul(cwave[ik], rt[ik][im]);
#endif
				}

				/* inverse Fourier transform */
				ifft(wave[im], cwavem);
			}

			/* left matrix */
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, im, c) shared(curr, it)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					i=iz+ix*nz;
					j=iz+ix*fnz;

					/* load source signature */
#ifdef SF_HAS_COMPLEX_H
					c=ww[it]*rr[i];
#else
					c=sf_crmul(ww[it], rr[i]);
#endif

					for(im=0; im<m; im++){
#ifdef SF_HAS_COMPLEX_H
						c += lt[im][i]*wave[im][j];
#else
						c = sf_cadd(c, sf_cmul(lt[im][i], wave[im][j]));
#endif
					}
					curr[j]=c;
				}
			}

			/* extract recorded data */
#ifdef _OPENMP
#pragma omp parallel for private(ir, i)
#endif
			for(ir=0; ir<nr; ir++){
				i=(ir*dr_v+r0_v)*fnz+rz;
				pp[ir][it]=curr[i];
			}

			if(!scomp){
				if(is==ns/2 && it%50==0){
					for(ix=0; ix<rnx; ix++){
						sf_complexwrite(curr+(ix+nb)*fnz+nb, rnz, Fwav);
					}
				}

				if(it%scalet==0 && function==3){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i)
#endif
					for(ix=0; ix<rnx; ix++){
						for(iz=0; iz<rnz; iz++){
							i=(ix+nb)*fnz+iz+nb;
							wfl[wfit][ix][iz]=curr[i];
						}
					}
					wfit++;
				}
			} //!scomp
		} // end of it 
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		/* compensated source wavefield */
		if(scomp){
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
			for(iz=0; iz<fnzx; iz++)
				curr[iz]=sf_cmplx(0.,0.);

			wfit=0;
			for(it=0; it<nt; it++){
				if(verb) sf_warning("Forward propagation it=%d/%d", it+1, nt);

				/* forward Fourier transform */
				fft(curr, cwave);

				/* right matrix */
				for(im=0; im<m2; im++){
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif 
					for(ik=0; ik<nk; ik++){
#ifdef SF_HAS_COMPLEX_H
						cwavem[ik]=cwave[ik]*rtb[ik][im];
#else
						cwavem[ik]=sf_cmul(cwave[ik], rtb[ik][im]);
#endif
					}

					/* inverse Fourier transform */
					ifft(wave2[im], cwavem);
				}

				/* left matrix */
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, im, c) shared(curr, it)
#endif
				for(ix=0; ix<nx; ix++){
					for(iz=0; iz<nz; iz++){
						i=iz+ix*nz;
						j=iz+ix*fnz;

						/* load source signature */
#ifdef SF_HAS_COMPLEX_H
						c=ww[it]*rr[i];
#else
						c=sf_crmul(ww[it], rr[i]);
#endif

						for(im=0; im<m2; im++){
#ifdef SF_HAS_COMPLEX_H
							c += ltb[im][i]*wave2[im][j];
#else
							c = sf_cadd(c, sf_cmul(ltb[im][i], wave2[im][j]));
#endif
						}
						curr[j]=c;
					}
				}

				if(is==ns/2 && it%50==0){
					for(ix=0; ix<rnx; ix++){
						sf_complexwrite(curr+(ix+nb)*fnz+nb, rnz, Fwav);
					}
				}

				if(it%scalet==0 && function==3){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i)
#endif
					for(ix=0; ix<rnx; ix++){
						for(iz=0; iz<rnz; iz++){
							i=(ix+nb)*fnz+iz+nb;
							wfl[wfit][ix][iz]=curr[i];
						}
					}
					wfit++;
				}
			} // end of it 
		} //scomp

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		/* check wfnt */
		if(function==3){ 
			if(wfit != wfnt) sf_error("At this moment, wfit should be equal to wfnt");
		}

		if(function==1){
			/* write data into disk */
			fseeko(swap, is*nr*nt*sizeof(float complex), SEEK_SET);
			fwrite(pp[0], sizeof(float complex)*nt, nr, swap);
			fclose(swap);
		}else{
			/* read data */
			sf_seek(Fdat, is*nr*nt*sizeof(float complex), SEEK_SET);
			sf_complexread(dd[0], nr*nt, Fdat);

			/* misfit calculation */
			for(ir=0; ir<nr; ir++){
				for(it=0; it<nt; it++){
					pp[ir][it] = dd[ir][it]-pp[ir][it];
					fcost += 0.5*cabsf(pp[ir][it])*cabsf(pp[ir][it]);
				}
			}

			if(is==ns/2) sf_complexwrite(pp[0], nr*nt, Fres);
		}

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		
		if(function==3){
			
			/* backward propagation */
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
			for(iz=0; iz<fnzx; iz++){
				curr[iz]=sf_cmplx(0.,0.);
			}
			
			wfit=wfnt-1;
			for(it=nt-1; it>=0; it--){ // it loop
				if(verb) sf_warning("Backward propagation it=%d/%d",it+1, nt);
            
				/* left matrix */
				for(im=0; im<m2; im++){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j) shared(curr, cwave)
#endif
					for(ix=0; ix<fnx; ix++){
						for(iz=0; iz<fnz; iz++){
							i=iz+ix*nz;
							j=iz+ix*fnz;

							if(ix<nx && iz<nz){
#ifdef SF_HAS_COMPLEX_H
								cwave[j]=conjf(ltb[im][i])*curr[j];
#else
								cwave[j]=sf_cmul(conjf(ltb[im][i]), curr[j]);
#endif
							}else{
								cwave[j]=sf_cmplx(0.,0.);
							}
						}
					}
					
					/* forward Fourier transform */
					fft(cwave, wave2[im]);
				}

#ifdef _OPENMP
#pragma omp parallel for private(ik, c, im)
#endif
				/* right matrix */
				for(ik=0; ik<nk; ik++){
					c=sf_cmplx(0., 0.);
					
					for(im=0; im<m2; im++){
#ifdef SF_HAS_COMPLEX_H
					c += wave2[im][ik]*conjf(rtb[ik][im]);
#else
					c = sf_cadd(c, sf_cmul(wave2[im][ik],conjf(rtb[ik][im])));
#endif
					}
					cwavem[ik]=c;
				}

				/* inverse Fourier transformm */
				ifft(curr, cwavem);

				/* load data residual */
#ifdef _OPENMP
#pragma omp parallel for private(ir,i)
#endif
				for(ir=0; ir<nr; ir++){
					i=(ir*dr_v+r0_v)*fnz+rz;
#ifdef SF_HAS_COMPLEX_H
					curr[i] += pp[ir][it];
#else
					curr[i]=sf_cadd(curr[i], pp[ir][it]);
#endif
				}

				if(is==ns/2 && it%50==0){
					for(ix=0; ix<rnx; ix++){
						sf_complexwrite(curr+(ix+nb)*fnz+nb, rnz, Fwav2);
					}
				}
				
				/* adjoint test */
				if(adjtest && is==ns/2){

					/* forward */
					for(ix=0; ix<nx; ix++){
						for(iz=0; iz<nz; iz++){
							i=iz+ix*nz;
							j=iz+ix*fnz;
							
							/* weighted source */
#ifdef SF_HAS_COMPLEX_H
							c=ww[it]*rr[i];
#else
							c=sf_crmul(ww[it], rr[i]);
#endif

							/* (f^*)*p */ 
#ifdef SF_HAS_COMPLEX_H
							forw += conjf(c)*curr[j];
#else
							forw=sf_cadd(forw, sf_cmul(conjf(c), curr[j]));
#endif
						}
					}

					/* adjoint */
					if(wfnt != nt) sf_error("need wfnt=nt for adjointness test");

					for(ir=0; ir<nr; ir++){
						i=ir*dr_v+r0_v-nb;
#ifdef SF_HAS_COMPLEX_H
						adjt += conjf(wfl[it][i][rz-nb])*pp[ir][it];
#else
						adjt = sf_cadd(adjt, sf_cmul(conjf(wfl[it][i][rz-nb]), pp[ir][it]));
#endif
					}
				}


				/* calculate gradient */
				if(it%scalet==0){
						fft(curr,cfrac);
						for(iz=0; iz<fnzx; iz++) {
#ifdef SF_HAS_COMPLEX_H
							cfrac1[iz]=cfrac[iz]*kk1[iz];
							cfrac2[iz]=cfrac[iz]*kk2[iz];
#else
							cfrac1[iz]=sf_crmul(cfrac[iz], kk1[iz]);
							cfrac2[iz]=sf_crmul(cfrac[iz], kk2[iz]);
#endif
						}
						ifft(cwave,cfrac1);
						ifft(cwavem,cfrac2);
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, q, v, a1, a2, wf)
#endif
						for(ix=0; ix<rnx; ix++){
							for(iz=0; iz<rnz; iz++){
								i  = ix*rnz+iz;
								j  = (ix+nb)*fnz+iz+nb;
								q  = qq[i];
								v  = vv[i];
								
								a1 = -(q+0.5)*powf(v, 2.*q)*powf(w0, -2.*q)*sinf(SF_PI*q)*powf(cosf(SF_PI*q/2.),2);
								a2 = (q+1.)*powf(v, q)*powf(w0, -q)*sqrtf(cosf(SF_PI*q))*cosf(SF_PI*q/2.);
								wf  = wfl[wfit][ix][iz];
								grad[i] += a1*crealf( conjf(wf) * cwave  [j] ) + 
									       a2*crealf( conjf(wf ) * sf_cmplx(0,1.)*cwavem [j] ); 
							}
						} 
					}
					wfit--;
			} // end of it
			if(wfit != -1) sf_error("Check program! The final wfit should be -1");

			if(adjtest && is==ns/2) sf_warning(" L[m]*d=(%g,%g) \n L'[d]*m=(%g,%g)", crealf(forw), cimagf(forw), crealf(adjt), cimagf(adjt));
		} // function==3

	} // end of shot
	MPI_Barrier(comm);

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	if(function==1){
		if(cpuid==0){
			swap=fopen("tempswap.bin","r");
			for(is=0; is<ns; is++){
				fseeko(swap, is*nr*nt*sizeof(float complex), SEEK_SET);
				fread(pp[0], sizeof(float complex), nr*nt, swap);
				sf_complexwrite(pp[0], nr*nt, Fdat);
			}
			fclose(swap);
			remove("tempswap.bin");
		}

		sf_fileclose(Fdat);
	}else{
		/* reduce misfit */
		if(cpuid==0){
			sendbuf=MPI_IN_PLACE;
			recvbuf=&fcost;
		}else{
			sendbuf=&fcost;
			recvbuf=NULL;
		}
		MPI_Reduce(sendbuf, recvbuf, 1, MPI_FLOAT, MPI_SUM, 0, comm);

		if(cpuid==0) sf_floatwrite(&fcost, 1, Fmisfit);
		sf_fileclose(Fmisfit);

		free(*dd); free(dd);
	}

	if(function==3){
		/* reduce gradient */
		if(cpuid==0){
			sendbuf=MPI_IN_PLACE;
			recvbuf=grad;
		}else{
			sendbuf=grad;
			recvbuf=NULL;
		}
		MPI_Reduce(sendbuf, recvbuf, rnzx, MPI_FLOAT, MPI_SUM, 0, comm);
		
		if(cpuid==0) sf_floatwrite(grad, rnzx, Fgrad);
		sf_fileclose(Fgrad);

		free(vv); free(qq); free(grad);
		free(**wfl); free(*wfl); free(wfl);
		free(kk1); free(kk2);
		free(cfrac); free(cfrac1); free(cfrac2);
		free(*ltb); free(ltb);
		free(*rtb); free(rtb);
		free(*wave2); free(wave2);
	}
	
	free(ww); free(rr); 
	free(*pp); free(pp);
	free(*lt); free(lt);
	free(*rt); free(rt);
		
	fft_finalize();
	free(*wave); free(wave);
	free(cwave); free(cwavem); free(curr);
	
	MPI_Finalize();
	exit(0);
}
