/* inverse Fourier transform of projection deviation operator */
/*
  Copyright (C) 2012 Tongji University (Jiubing Cheng)
 
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
#include "_cjb.h"

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

void kxkz2xz(float**xi, float **xo, int hnkx, int hnkz, int nkx, int nkz)
/*< kxkz2xz: inverse Fourier transform of operator from (kx,kz) to (x, z) domain>*/
{
       int ii, jj, i, j, nkxz;
       int subx, subz;
       sf_complex *xin, *xout;
#ifdef SF_HAS_FFTW
       fftwf_plan xpi;
#endif

       nkxz=nkx*nkz;

       xin=sf_complexalloc(nkxz);
       xout=sf_complexalloc(nkxz);

#ifdef SF_HAS_FFTW
       xpi=fftwf_plan_dft_2d(nkx,nkz,
			    (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);
#endif

	for(i=-hnkx;i<=hnkx;i++)
	for(j=-hnkz;j<=hnkz;j++)
	{
                ii=i+hnkx;
                jj=j+hnkz;

                if(i<=0)
                    subx=ii+hnkx;
                else
                    subx=ii-hnkx-1;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xin[ii*nkz+jj]=sf_cmplx(xi[subx][subz],0.);
	}

#ifdef SF_HAS_FFTW
	fftwf_execute(xpi);
#endif

        for(i=-hnkx;i<=hnkx;i++)
        for(j=-hnkz;j<=hnkz;j++)
        {
                ii=i+hnkx;
                jj=j+hnkz;

                if(i<=0)
                    subx=ii+hnkx;
                else
                    subx=ii-hnkx-1;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xo[subx][subz]=crealf(xout[ii*nkz+jj])/nkxz;
	}

#ifdef SF_HAS_FFTW
        fftwf_destroy_plan(xpi);
#endif

        free(xin);
        free(xout);
}

void ikxkz2xz(float**xi, float **xo, int hnkx, int hnkz, int nkx, int nkz)
/*< ikxkz2xz: inverse Fourier transform of operator (*i) from (kx,kz) to (x, z) domain>*/
{
       int ii, jj, i, j, nkxz;
       int subx, subz;
       sf_complex *xin, *xout;
#ifdef SF_HAS_FFTW
       fftwf_plan xpi;
#endif

       nkxz=nkx*nkz;

       xin=sf_complexalloc(nkxz);
       xout=sf_complexalloc(nkxz);

#ifdef SF_HAS_FFTW
       xpi=fftwf_plan_dft_2d(nkx,nkz,
			     (fftwf_complex *) xin, (fftwf_complex *) xout,
			     FFTW_BACKWARD,FFTW_ESTIMATE);
#endif

       for(i=-hnkx;i<=hnkx;i++)
       for(j=-hnkz;j<=hnkz;j++)
	{
                ii=i+hnkx;
                jj=j+hnkz;

                if(i<=0)
                    subx=ii+hnkx;
                else
                    subx=ii-hnkx-1;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xin[ii*nkz+jj]=sf_cmplx(0.,xi[subx][subz]);
	}

#ifdef SF_HAS_FFTW
	fftwf_execute(xpi);
#endif

        for(i=-hnkx;i<=hnkx;i++)
        for(j=-hnkz;j<=hnkz;j++)
        {
                ii=i+hnkx;
                jj=j+hnkz;

                if(i<=0)
                    subx=ii+hnkx;
                else
                    subx=ii-hnkx-1;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xo[subx][subz]=crealf(xout[ii*nkz+jj])/(nkxz);
	}

#ifdef SF_HAS_FFTW
        fftwf_destroy_plan(xpi);
#endif

        free(xin);
        free(xout);
}
