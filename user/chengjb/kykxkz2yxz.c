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
       int ii, jj, i, j, iinkz, nkxz;
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
        {
            ii=i+hnkx;
            if(i<=0)
                subx=ii+hnkx;
            else
                subx=ii-hnkx-1;

            iinkz=ii*nkz;

	    for(j=-hnkz;j<=hnkz;j++)
	    {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xin[iinkz+jj]=sf_cmplx(xi[subx][subz],0.);
            }
	}

#ifdef SF_HAS_FFTW
	fftwf_execute(xpi);
#endif

        for(i=-hnkx;i<=hnkx;i++)
        {
             ii=i+hnkx;
             if(i<=0)
                subx=ii+hnkx;
             else
                subx=ii-hnkx-1;

             iinkz=ii*nkz;

             for(j=-hnkz;j<=hnkz;j++)
             {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xo[subx][subz]=crealf(xout[iinkz+jj])/nkxz;
            }
	}

#ifdef SF_HAS_FFTW
        fftwf_destroy_plan(xpi);
#endif

        free(xin);
        free(xout);
}

void kykxkz2yxz(float***xi, float ***xo, int hnky, int hnkx, int hnkz, int nky, int nkx, int nkz)
/*< kykxkz2yxz: inverse Fourier transform of operator from (ky,kx,kz) to (y, x, z) domain>*/
{
       int ii, jj, kk, i, j, k, nkxz, nkxyz;
       int kknkxz, iinkz, kij;
       int suby, subx, subz;

       sf_complex *xin, *xout;
#ifdef SF_HAS_FFTW
       fftwf_plan xpi;
#endif

       nkxz=nkx*nkz;
       nkxyz=nky*nkx*nkz;

       xin=sf_complexalloc(nkxyz);
       xout=sf_complexalloc(nkxyz);

#ifdef SF_HAS_FFTW
       xpi=fftwf_plan_dft_3d(nky, nkx, nkz,
			    (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);
#endif
      for(k=-hnky;k<=hnky;k++)
      {
        kk=k+hnky;
        if(k<=0)
            suby=kk+hnky;
        else
            suby=kk-hnky-1;

        kknkxz = kk*nkxz;

	for(i=-hnkx;i<=hnkx;i++)
        {
            ii=i+hnkx;
            if(i<=0)
               subx=ii+hnkx;
            else
               subx=ii-hnkx-1;

            iinkz=ii*nkz;
            kij=kknkxz+iinkz;

	    for(j=-hnkz;j<=hnkz;j++)
	    {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xin[kij+jj]=sf_cmplx(xi[suby][subx][subz],0.);
            }
         }
      }

#ifdef SF_HAS_FFTW
      fftwf_execute(xpi);
#endif
      for(k=-hnky;k<=hnky;k++)
      {
        kk=k+hnky;
        if(k<=0)
            suby=kk+hnky;
        else
            suby=kk-hnky-1;

        kknkxz = kk*nkxz;

        for(i=-hnkx;i<=hnkx;i++)
        {
            ii=i+hnkx;
            if(i<=0)
                subx=ii+hnkx;
            else
                subx=ii-hnkx-1;

            iinkz=ii*nkz;
            kij=kknkxz+iinkz;

            for(j=-hnkz;j<=hnkz;j++)
            {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xo[suby][subx][subz]=crealf(xout[kij+jj])/nkxyz;
            }
	 }
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
       int ii, jj, i, j, nkxz, iinkz;
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
       {
            ii=i+hnkx;
            if(i<=0)
                subx=ii+hnkx;
            else
                subx=ii-hnkx-1;

            iinkz=ii*nkz;

            for(j=-hnkz;j<=hnkz;j++)
	    {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xin[iinkz+jj]=sf_cmplx(0.,xi[subx][subz]);
            }
	}

#ifdef SF_HAS_FFTW
	fftwf_execute(xpi);
#endif

        for(i=-hnkx;i<=hnkx;i++)
        {
             ii=i+hnkx;
             if(i<=0)
                subx=ii+hnkx;
             else
                subx=ii-hnkx-1;

             iinkz=ii*nkz;

             for(j=-hnkz;j<=hnkz;j++)
             {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xo[subx][subz]=crealf(xout[iinkz+jj])/(nkxz);
             }
	}

#ifdef SF_HAS_FFTW
        fftwf_destroy_plan(xpi);
#endif

        free(xin);
        free(xout);
}

void ikykxkz2yxz(float***xi, float ***xo, int hnky, int hnkx, int hnkz, int nky, int nkx, int nkz)
/*< ikykxkz2yxz: inverse Fourier transform of operator (*i) from (ky, kx,kz) to (y, x, z) domain>*/
{
       int ii, jj, kk, i, j, k, nkxz, nkxyz;
       int kknkxz, iinkz, kij;
       int suby, subx, subz;
       sf_complex *xin, *xout;
#ifdef SF_HAS_FFTW
       fftwf_plan xpi;
#endif

       nkxz=nkx*nkz;
       nkxyz=nky*nkx*nkz;

       xin=sf_complexalloc(nkxyz);
       xout=sf_complexalloc(nkxyz);

#ifdef SF_HAS_FFTW
       xpi=fftwf_plan_dft_3d(nky, nkx, nkz,
			     (fftwf_complex *) xin, (fftwf_complex *) xout,
			     FFTW_BACKWARD,FFTW_ESTIMATE);
#endif
     for(k=-hnky;k<=hnky;k++)
     {
       kk=k+hnky;
       if(k<=0)
           suby=kk+hnky;
       else
           suby=kk-hnky-1;

        kknkxz = kk*nkxz;

       for(i=-hnkx;i<=hnkx;i++)
       {
            ii=i+hnkx;
            if(i<=0)
                subx=ii+hnkx;
            else
                subx=ii-hnkx-1;

            iinkz=ii*nkz;
            kij=kknkxz+iinkz;

            for(j=-hnkz;j<=hnkz;j++)
	    {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xin[kij+jj]=sf_cmplx(0.,xi[suby][subx][subz]);
            }
	}
     }

#ifdef SF_HAS_FFTW
     fftwf_execute(xpi);
#endif
     for(k=-hnky;k<=hnky;k++)
     {
        kk=k+hnky;
        if(k<=0)
           suby=kk+hnky;
        else
           suby=kk-hnky-1;

        kknkxz = kk*nkxz;

        for(i=-hnkx;i<=hnkx;i++)
        {
             ii=i+hnkx;
             if(i<=0)
                subx=ii+hnkx;
             else
                subx=ii-hnkx-1;

             iinkz=ii*nkz;
             kij=kknkxz+iinkz;

             for(j=-hnkz;j<=hnkz;j++)
             {
                jj=j+hnkz;

                if(j<=0)
                    subz=jj+hnkz;
                else
                    subz=jj-hnkz-1;

		xo[suby][subx][subz]=crealf(xout[kij+jj])/(nkxyz);
             }
	 }
       }

#ifdef SF_HAS_FFTW
       fftwf_destroy_plan(xpi);
#endif

       free(xin);
       free(xout);
}
