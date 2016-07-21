/* Complex-valued compensation (between two wavefields) */
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "cfft2w.h"
#include "compensate.h"

int main(int argc, char* argv[])
{
    bool verb,cmplx;
    int i,j,ix,iz,iter,niter; /* index variables */
    int nt,nz,nx,nzx,nz2,nx2,nzx2,nk,pad1,nth;
    float perc,vmax,eps;

    /* I/O arrays*/
    sf_complex *num0,*num,*num2,*den0,*den,*wht,*iwht,*fnum,*fden,*fden2,*fwht;
    /* I/O files */
    sf_file Fnum,Fden,Fres;
    /* cube axes */
    sf_axis at,az,ax;

    sf_init(argc,argv);
    if(!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    if(!sf_getbool("cmplx",&cmplx)) cmplx=true; /* use complex i/o */
    if(!sf_getint("niter",&niter)) niter=1; /* number of iterations */
    if(!sf_getfloat("perc",&perc)) perc=0.1; /* precentage (of max) for protection when dividing */
    perc /= 100.;

    /* setup I/O files */
    Fnum = sf_input("in" );
    Fden = sf_input("den");
    Fres = sf_output("out");

    if (SF_COMPLEX != sf_gettype(Fnum)) sf_error("Need complex input");
    if (SF_COMPLEX != sf_gettype(Fden)) sf_error("Need complex den");

    sf_settype(Fres,SF_COMPLEX);

    /* Read/Write axes */
    az = sf_iaxa(Fnum,1); nz = sf_n(az); 
    ax = sf_iaxa(Fnum,2); nx = sf_n(ax); 
    at = sf_iaxa(Fnum,3); nt = sf_n(at); 

    sf_oaxa(Fres,az,1); 
    sf_oaxa(Fres,ax,2); 
    sf_oaxa(Fres,at,3);
    
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

#ifdef _OPENMP
#pragma omp parallel
    {
      nth = omp_get_num_threads();
    }
    if (verb) sf_warning(">>>> Using %d threads <<<<<", nth);
#endif

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);

    nzx = nz*nx;
    nzx2 = nz2*nx2;

    /* Input wavefields */
    num0 = sf_complexalloc(nzx );
    num  = sf_complexalloc(nzx2);
    num2 = sf_complexalloc(nzx2);
    den0 = sf_complexalloc(nzx );
    den  = sf_complexalloc(nzx2);
    wht  = sf_complexalloc(nzx2);
    iwht = sf_complexalloc(nzx2);
    fnum = sf_complexalloc(nk  );
    fden = sf_complexalloc(nk  );
    fden2= sf_complexalloc(nk  );
    fwht = sf_complexalloc(nk  );

    sf_complexread(num0,nzx,Fnum);
    sf_complexread(den0,nzx,Fden);

    sf_fileclose(Fnum);
    sf_fileclose(Fden);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < nzx2; i++) {
        wht[i] = sf_cmplx(1.f,0.f);
        num[i] = sf_cmplx(0.f,0.f);
        den[i] = sf_cmplx(0.f,0.f);
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < nk; i++) fwht[i] = sf_cmplx(1.f,0.f);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,iz,i,j)
#endif
    for (ix = 0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            i = iz+ix*nz;  /* original grid */
            j = iz+ix*nz2; /* padded grid */
            num[j] = num0[i];
            den[j] = den0[i];
        }
    }

    /***************************************/
    /* start the work */

    cfft2(den,fden);

    for (iter=0; iter<niter; iter++) {
	if(verb) sf_warning("iter=%d/%d;",iter,niter);

        /* apply Fourier scaling to denominator */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
        for (i = 0; i < nk; i++) {
#ifdef SF_HAS_COMPLEX_H
            fden2[i] = fden[i]*fwht[i];
#else
            fden2[i] = sf_cmul(fden[i],fwht[i]);
#endif
        }

        icfft2(den,fden2);

        /* stable division in the space domain: num/den */
        vmax = find_cmax(nzx, den);
        eps = (vmax==0) ? SF_EPS : vmax*vmax*perc;
        if (verb) sf_warning("Space domain: vmax=%f, eps=%f",vmax,eps);
        stable_cdiv(nzx, eps, num, den, wht);

        /* apply inverse weight to numerator */
        cinvert(nzx, wht, iwht);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
        for (i = 0; i < nzx; i++) {
#ifdef SF_HAS_COMPLEX_H
            num2[i] = num[i]*iwht[i];
#else
            num2[i] = sf_cmul(num[i],iwht[i]);
#endif
        }

        /* forward FFT */
        cfft2(num2,fnum);
        //cfft2(den,fden);

        vmax = find_cmax(nk, fden);
        eps = (vmax==0) ? SF_EPS : vmax*vmax*perc;
        if (verb) sf_warning("Fourier domain: vmax=%f, eps=%f",vmax,eps);
        stable_cdiv(nk, eps, fnum, fden, fwht);
    }
    if(verb) sf_warning("."); 

    /* apply scaling */
    cfft2(num,fnum);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < nk; i++) {
#ifdef SF_HAS_COMPLEX_H
        fnum[i] = fnum[i]*fwht[i];
#else
        fnum[i] = sf_cmul(fnum[i],fwht[i]);
#endif
    }

    icfft2(num2,fnum);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < nzx; i++) {
#ifdef SF_HAS_COMPLEX_H
        num2[i] = num2[i]*wht[i];
#else
        num2[i] = sf_cmul(num2[i],wht[i]);
#endif
    }

    /* output result */
    sf_complexwrite(num2,nzx,Fres);

    cfft2_finalize();
    exit (0);
}
