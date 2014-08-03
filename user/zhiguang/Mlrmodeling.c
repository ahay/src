/* 2-D FFT-based modeling */
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

#include "fft2.h"

int main(int argc, char* argv[])
{
	bool wantwf, cmplx;
	int ix, iz, it, i, j, in, ik;
	int nt, nx, nz, nx2, nz2, nzx, nzx2, nk, pad1, n2, m2, snpint, wfnt;
	int n0, nth;
	float dt, t0, dz, z0, old, c, x0, dx, wfdt;
	/* boundary condition variables */
	int nb, rnz, rnx;
	float *bcx, *bcz, coef;
	
	float *ww, *rr, **dd, **lt, **rt;
	float *curr, *prev, **wave;
	sf_complex *cwave, *cwavem;

	sf_axis at, ax, az;
	sf_file Frcd, Fsrc, Frr, Fleft, Fright, Fwf;

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
	if(!sf_getint("pad1", &pad1)) pad1=1;
	if(!sf_getint("nb", &nb)) nb=40; /* padded boundary width */
	if(!sf_getfloat("coef", &coef)) coef=0.01; /* decaying parameter */
	if(!sf_getint("n0", &n0))   sf_error("Need n0=  !");
	if(!sf_getint("snapint", &snpint)) snpint=10; /* Just for output wavefield */

	Fsrc=sf_input("in");
	Frcd=sf_output("out");
	Frr=sf_input("ref");

	sf_warning("Just a mark!");
	at=sf_iaxa(Fsrc, 1); nt=sf_n(at); dt=sf_d(at); t0=sf_o(at);
	ax=sf_iaxa(Frr, 2); nx=sf_n(ax); dx=sf_d(ax); x0=sf_o(ax);
	az=sf_iaxa(Frr, 1); nz=sf_n(az); dz=sf_d(az); z0=sf_o(az);
	rnx=nx-2*nb; sf_setn(ax, rnx); sf_seto(ax, x0+dx*nb);
	rnz=nz-2*nb; sf_setn(az, rnz); sf_seto(az, z0+dz*nb);
	sf_oaxa(Frcd, at, 1);
	sf_oaxa(Frcd, ax, 2);

	ww=sf_floatalloc(nt);
	sf_floatread(ww, nt, Fsrc);
	dd=sf_floatalloc2(nt, rnx); 
	
	Fleft=sf_input("left");
	Fright=sf_input("right");

	sf_warning("nx=%d, nz=%d, x0=%g, z0=%g, rnx=%d, rnz=%d",nx,nz,x0,z0,rnx,rnz);
	nk=fft2_init(cmplx, pad1, nz, nx, &nz2, &nx2);
	nzx=nz*nx;
	nzx2=nz2*nx2;
	if(!sf_histint(Fleft, "n1", &n2) || n2!=nzx) sf_error("Need n1=%d in left", nzx);
	if(!sf_histint(Fleft, "n2", &n2)) sf_error("No n2= in left");
	if(!sf_histint(Fright, "n1", &m2) || m2!= n2) sf_error("Need n1=%d in right", n2);
	if(!sf_histint(Fright, "n2", &m2) || m2!= nk) sf_error("Need n2=%d in right", nk);

	rr=sf_floatalloc(nzx);
	sf_floatread(rr, nzx, Frr);

	lt=sf_floatalloc2(nzx, n2);
	rt=sf_floatalloc2(n2, nk);
	sf_floatread(lt[0], nzx*n2, Fleft);
	sf_floatread(rt[0], n2*nk, Fright);

	sf_warning("Just a mark!");
	curr=sf_floatalloc(nzx2);
	prev=sf_floatalloc(nzx);
	cwave=sf_complexalloc(nk);
	cwavem=sf_complexalloc(nk);
	wave=sf_floatalloc2(nzx2, n2);

	if(wantwf){
		wfnt=(nt-1)/snpint+1;
		wfdt=dt*snpint;

		Fwf=sf_output("wf");
		sf_oaxa(Fwf, az, 1);
		sf_oaxa(Fwf, ax, 2);
		sf_putint(Fwf, "n3", wfnt);
		sf_putfloat(Fwf, "d3", wfdt);
		sf_putfloat(Fwf, "o3", t0);
	}

	/* Calculate absorption coefficients */
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

	ifft2_allocate(cwavem);

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
	
	for(it=0; it<nt; it++){
		sf_warning("Source(VP) it=%d/%d;", it+1,nt);

		for(ix=0; ix<rnx; ix++)
		dd[ix][it]=curr[(ix+nb)*nz2+n0+nb];

		fft2(curr, cwave);

		for(in=0; in<n2; in++){
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
			for(ik=0; ik<nk; ik++){
#ifdef SF_HAS_COMPLEX_H
				cwavem[ik]=cwave[ik]*rt[ik][in];
#else
				cwavem[ik]=sf_crmul(cwave[ik],rt[ik][in]);
#endif
			}
			ifft2(wave[in],cwavem);
		}

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz, i, j, old, c, in) shared(curr, prev, rr, lt, wave)
#endif
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				i=ix*nz+iz; //original grid
				j=ix*nz2+iz; //padded grid
				old=c=curr[j];
				c += c + ww[it]*rr[i] - prev[i];
				prev[i] = old*bcx[ix]*bcz[iz];

				for(in=0; in<n2; in++){
					c+=lt[in][i]*wave[in][j];
				}
				curr[j] = c*bcx[ix]*bcz[iz];
			}
		}

		if(wantwf && it%snpint==0){
			for(ix=0; ix<rnx; ix++){
				sf_floatwrite(curr+(ix+nb)*nz2+nb, rnz, Fwf);
			}
		}
	} //end of it

	sf_floatwrite(dd[0], rnx*nt, Frcd);

	exit(0);
}
