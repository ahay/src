/* 2-D FFT-based RTM for converted wave */
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
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fft2.h"

int main(int argc, char* argv[])
{
	bool wantwf, cmplx, correct;
	int ix, iz, it, is, wfit, i, j, in, ik;
	int nt, nt2, wfnt, snpint, nx, nz, nx2, nz2, nvx, nzx, nzx2, nk, pad1, n2, m2, tmpint;
	int n0, nth, multiple, snpint2, snapshot, ns, move;
	float dt, wfdt, t0, old, c, ds, s0, tmpfloat, dx, x0;
	/* boundary condition variables */
	int nb, rnz, rnx, nc;
	float *bcx, *bcz, coef, frequency;
	char *path1, *path2, number[5], *leftf, *leftb, *rightf, *rightb;
	
	float *ww, *rr, **vmm, ***wf2; /* No change with different shot number */
	float *curr, *prev, **dd, **mm, ***wf1; /* Change with shot number */
	float **ltf, **rtf, **ltb, **rtb, **wave, **waveb; /* Even with different sizes */
	sf_complex *cwave, *cwavem;

	sf_axis at, ax, az, as;
	sf_file Fvv, Frcd, Fimg, Fsrc, Frr, Fleftf, Frightf, Fleftb, Frightb, Fwf1, Fwf2;

	sf_init(argc, argv);

#ifdef _OPENMP
#pragma omp parallel
	{
		nth = omp_get_num_threads();
	}
	sf_warning(">>> Using %d threads <<<", nth);
#endif

	if(!sf_getbool("wantwf", &wantwf)) wantwf=false;
	if(!sf_getbool("cmplx", &cmplx)) cmplx=false;
	if(!sf_getbool("correct", &correct)) correct=false;
	if(!sf_getint("pad1", &pad1)) pad1=1;
	if(!sf_getint("nb", &nb)) nb=40; /* padded boundary width */
	if(!sf_getfloat("coef", &coef)) coef=0.01; /* decaying parameter */
	if(!sf_getint("multiple", &multiple)) multiple=2;
	if(!sf_getfloat("frequency", &frequency)) frequency=10.0;
	if(!sf_getint("snapshot", &snapshot)) snapshot=1; /* The default is the first shot */
	if(!sf_getint("n0", &n0))   sf_error("Need n0=  !"); /* n0 start from the padded model */
	if(!sf_getint("snapint", &snpint)) snpint=10;

	Fvv=sf_input("in");
	Fimg=sf_output("out");
	Frcd=sf_input("data");
	Fsrc=sf_input("src");
	Frr=sf_input("ref");

	at=sf_iaxa(Fsrc, 1); nt=sf_n(at); dt=sf_d(at); t0=sf_o(at);
	ax=sf_iaxa(Frcd, 2); rnx=sf_n(ax); dx=sf_d(ax); x0=sf_o(ax); nx=rnx+2*nb;
	az=sf_iaxa(Fvv, 1); rnz=sf_n(az); nz=rnz+2*nb;
	as=sf_iaxa(Frcd, 3); ns=sf_n(as); ds=sf_d(as); s0=sf_o(as);
	move=(int)(ds/dx);

	/* check input and output file dimensions */
	nt2=(nt-1)/multiple+1;
	if(!sf_histint(Frcd, "n1", &tmpint) || tmpint != nt2) sf_error("n1 in data should be %d",nt2);
	/* check time axis */
	if(!sf_histint(Frr, "n1", &tmpint) || tmpint != nz) sf_error("n1 in ref should be %d", nz);
	if(!sf_histint(Frr, "n2", &tmpint) || tmpint != nx) sf_error("n2 in ref should be %d", nx);
	/* check padding dimensions */
	if(!sf_histint(Fvv, "n2", &nvx) || nvx != (int)(ds/dx*(ns-1)+rnx)) sf_error("Wrong dimensions!");
	if(!sf_histfloat(Fvv, "o2", &tmpfloat) || tmpfloat != s0+x0) sf_error("Wrong dimensions!");
	/* check Fvv */
	sf_warning("ns=%d, nvx=%d, nx=%d, nz=%d, rnx=%d, rnz=%d, move=%d",ns,nvx,nx,nz,rnx,rnz,move);

	nk=fft2_init(cmplx, pad1, nz, nx, &nz2, &nx2);
	nzx=nz*nx;
	nzx2=nz2*nx2;

	/* read wavelet and source map */
	rr=sf_floatalloc(nzx);
	sf_floatread(rr, nzx, Frr);
	ww=sf_floatalloc(nt);
	sf_floatread(ww, nt, Fsrc);

	/* storage allocation */
	dd=sf_floatalloc2(nt2, rnx);
	vmm=sf_floatalloc2(rnz, nvx);
	mm=sf_floatalloc2(rnz, rnx);
	for(ix=0; ix<nvx; ix++)
		for(iz=0; iz<rnz; iz++)
			vmm[ix][iz]=0.0;
	
	curr=sf_floatalloc(nzx2);
	prev=sf_floatalloc(nzx);
	cwave=sf_complexalloc(nk);
	cwavem=sf_complexalloc(nk);

	/* wavefield parameters */
	snpint2=snpint/multiple;
	wfnt=(nt-1)/snpint+1;
	wfdt=dt*snpint;
	wf1=sf_floatalloc3(rnz, rnx, wfnt);
	if(wantwf){
		Fwf1=sf_output("wf1");
		Fwf2=sf_output("wf2");
		wf2=sf_floatalloc3(rnz, rnx, wfnt);

		sf_oaxa(Fwf1, az, 1);
		sf_oaxa(Fwf1, ax, 2);
		sf_putint(Fwf1, "n3", wfnt);
		sf_putfloat(Fwf1, "d3", wfdt);
		sf_putfloat(Fwf1, "o3", t0);

		sf_oaxa(Fwf2, az, 1);
		sf_oaxa(Fwf2, ax, 2);
		sf_putint(Fwf2, "n3", wfnt);
		sf_putfloat(Fwf2, "d3", -wfdt);
		sf_putfloat(Fwf2, "o3", (wfnt-1)*wfdt);
	}

	if(correct){
		nc=(int)(1.0/frequency/wfdt+0.5);
	}else{
		nc=0;
	}
	sf_warning("Wavefield correction nc=%d", nc);

	/* calculate absorption coefficients */
	bcx=sf_floatalloc(nx);
	bcz=sf_floatalloc(nz);
	for(ix=0; ix<nx; ix++)
		bcx[ix]=1.0;
	for(iz=0; iz<nz; iz++)
		bcz[iz]=1.0;
	for(ix=0; ix<nb; ix++){
		bcx[ix]=exp(-coef*coef*(nb-1-ix)*(nb-1-ix));
		bcz[ix]=bcx[ix];
		bcx[nx-1-ix]=bcx[ix];
		bcz[nz-1-ix]=bcz[ix];
	}

	path1=sf_getstring("path1");
	path2=sf_getstring("path2");
	if(path1==NULL) path1="./mat/left";
	if(path2==NULL) path2="./mat/right";
	/* shot loop */
	for(is=0; is<ns; is++){
		/* construct the names of left and right matrices */
		leftf=sf_charalloc(strlen(path1));
		rightf=sf_charalloc(strlen(path2));
		leftb=sf_charalloc(strlen(path1));
		rightb=sf_charalloc(strlen(path2));
		strcpy(leftf,path1);
		strcpy(rightf,path2);
		strcpy(leftb,path1);
		strcpy(rightb,path2);
		sprintf(number,"%d",is+1);
		
		strcat(leftf,number);
		strcat(rightf,number);
		strcat(leftb,number);
		strcat(rightb,number);
		strcat(leftf,"f");
		strcat(rightf,"f");
		strcat(leftb,"b");
		strcat(rightb,"b");
		
		Fleftf=sf_input(leftf);
		Frightf=sf_input(rightf);
		Fleftb=sf_input(leftb);
		Frightb=sf_input(rightb);
		
		/* read dimensions */
		if(!sf_histint(Fleftf, "n1", &n2) || n2!=nzx) sf_error("Need n1=%d in leftf", nzx);
		if(!sf_histint(Fleftf, "n2", &n2)) sf_error("No n2= in leftf");
		if(!sf_histint(Frightf, "n1", &m2) || m2!= n2) sf_error("Need n1=%d in rightf", n2);
		if(!sf_histint(Frightf, "n2", &m2) || m2!= nk) sf_error("Need n2=%d in rightf", nk);
		if(!sf_histint(Fleftb, "n2", &m2)) sf_error("No n2= in leftb");
		if(!sf_histint(Frightb, "n1", &tmpint) || tmpint != m2) sf_error("Need n1=%d in rightb", m2);

		/* storage allocation for each shot migration */
		ltf=sf_floatalloc2(nzx, n2);
		rtf=sf_floatalloc2(n2, nk);
		ltb=sf_floatalloc2(nzx, m2);
		rtb=sf_floatalloc2(m2, nk);
		sf_floatread(ltf[0], nzx*n2, Fleftf);
		sf_floatread(rtf[0], n2*nk, Frightf);
		sf_floatread(ltb[0], nzx*m2, Fleftb);
		sf_floatread(rtb[0], m2*nk, Frightb);
		
		sf_warning("Just a mark!");
		wave=sf_floatalloc2(nzx2, n2);
		waveb=sf_floatalloc2(nzx2, m2);
		
		ifft2_allocate(cwavem);
		
		/* read data */
		sf_floatread(dd[0], nt2*rnx, Frcd);
	
		// Source wavefield propagation

#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
		for(iz=0; iz<nzx; iz++){
			prev[iz]=0.;
		}
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
		for(iz=0; iz<nzx2; iz++){
			curr[iz]=0.;
		}
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
		for(ix=0; ix<rnx; ix++){
			for(iz=0; iz<rnz; iz++){
				mm[ix][iz]=0.;
			}
		}
		
		wfit=0;
		for(it=0; it<nt; it++){
			sf_warning("Shot: %d; Source(VP) it=%d/%d;", is+1, it+1,nt);
			
			fft2(curr, cwave);
			
			for(in=0; in<n2; in++){
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
				for(ik=0; ik<nk; ik++){
#ifdef SF_HAS_COMPLEX_H
					cwavem[ik]=cwave[ik]*rtf[ik][in];
#else
					cwavem[ik]=sf_crmul(cwave[ik],rtf[ik][in]);
#endif
				}
				ifft2(wave[in],cwavem);
			}

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, old, c, in) shared(curr, prev, rr, ltf, wave)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					i=ix*nz+iz; //original grid
					j=ix*nz2+iz; //padded grid
					old=c=curr[j];
					c += c + ww[it]*rr[i] - prev[i];
					prev[i] = old*bcx[ix]*bcz[iz];
					
					for(in=0; in<n2; in++){
						c+=ltf[in][i]*wave[in][j];
					}
					curr[j] = c*bcx[ix]*bcz[iz];
				}
			}
			
			if(it%snpint==0){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, j)
#endif
				for(ix=0; ix<rnx; ix++){
					for(iz=0; iz<rnz; iz++){
						j=(ix+nb)*nz2+iz+nb;
						wf1[wfit][ix][iz]=curr[j];
					}
				}
				wfit++;
			}
		} //end of it
		
		// Receiver wavefield propagation
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
		for(iz=0; iz<nzx; iz++){
			prev[iz]=0.;
		}
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
		for(iz=0; iz<nzx2; iz++){
			curr[iz]=0.;
		}
		
		for(it=nt2-1; it>=0; it--){
			sf_warning("Shot: %d; Receiver(SV) it=%d/%d;", is+1, it+1, nt2);

#ifdef _OPENMP
#pragma omp parallel for private(ix, j) shared(n0, it)
#endif
			for(ix=0; ix<rnx; ix++){
				j=(ix+nb)*nz2+n0;
				curr[j]+=dd[ix][it];
			}
			
			fft2(curr, cwave);
			
			for(in=0; in<m2; in++){
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
				for(ik=0; ik<nk; ik++){
#ifdef SF_HAS_COMPLEX_H
					cwavem[ik]=cwave[ik]*rtb[ik][in];
#else
					cwavem[ik]=sf_crmul(cwave[ik],rtb[ik][in]);
#endif
				}
				
				ifft2(waveb[in], cwavem);
			}

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, old, c, in) shared(curr, prev, ltb, waveb)
#endif
			for(ix=0; ix<nx; ix++){
				for(iz=0; iz<nz; iz++){
					i=ix*nz+iz; //original grid
					j=ix*nz2+iz; //padded grid
					old=c=curr[j];
					c += c - prev[i];
					prev[i]=old*bcx[ix]*bcz[iz];
					
					for(in=0; in<m2; in++){
						c += ltb[in][i]*waveb[in][j];
					}
					
					curr[j]=c*bcx[ix]*bcz[iz];
				}
			}
			
			if(it%snpint2==0){
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, j)
#endif
				for(ix=0; ix<rnx; ix++){
					for(iz=0; iz<rnz; iz++){
						j=(ix+nb)*nz2+iz+nb;
						if(wfit+nc-1<wfnt) mm[ix][iz]+=curr[j]*wf1[wfit-1+nc][ix][iz];
						if(wantwf && is+1 == snapshot) wf2[wfnt-wfit][ix][iz]=curr[j];
					}
				}
				wfit--;
			}
		} //end of it
		
		if(is+1 == snapshot && wantwf){
			sf_floatwrite(wf1[0][0], wfnt*rnx*rnz, Fwf1);
			sf_floatwrite(wf2[0][0], wfnt*rnx*rnz, Fwf2);
			sf_fileclose(Fwf1);
			sf_fileclose(Fwf2);
		}
		
		for(ix=0; ix<rnx; ix++)
			for(iz=0; iz<rnz; iz++)
				vmm[ix+is*move][iz]+=mm[ix][iz];
		
		/* free memeory */
		free(*ltf); free(ltf);
		free(*rtf); free(rtf);
		free(*ltb); free(ltb);
		free(*rtb); free(rtb);
		free(*wave); free(wave);
		free(*waveb); free(waveb);
		free(leftf); free(rightf);
		free(leftb); free(rightb);
	}

	sf_floatwrite(vmm[0], nvx*rnz, Fimg);

	exit(0);
}
