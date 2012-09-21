/*************************************************************************
 * * inverse Fourier transform of projection deviation operator
 * *
 * *    Copyright: Tongji University (Jiubing Cheng)
 * *    2012.3.2
 * *************************************************************************/
//#include <rsf.h>
#include "_cjb.h"

#include <fftw3.h>

void kxkz2xz(float**xi, float **xo, int hnkx, int hnkz, int nkx, int nkz)
/*< kxkz2xz: inverse Fourier transform of operator from (kx,kz) to (x, z) domain>*/
{
       int ii, jj, i, j, nkxz;
       int subx, subz;

       nkxz=nkx*nkz;

       fftw_complex* xin=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nkxz);
       fftw_complex* xout=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nkxz);

       fftw_plan xpi=fftw_plan_dft_2d(nkx,nkz,xin,xout,FFTW_BACKWARD,FFTW_ESTIMATE);

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

		xin[ii*nkz+jj]=xi[subx][subz];
	}

	fftw_execute(xpi);

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

		xo[subx][subz]=xout[ii*nkz+jj]/nkxz;
	}

        fftw_destroy_plan(xpi);

        fftw_free(xin);
        fftw_free(xout);
}

void ikxkz2xz(float**xi, float **xo, int hnkx, int hnkz, int nkx, int nkz)
/*< ikxkz2xz: inverse Fourier transform of operator (*i) from (kx,kz) to (x, z) domain>*/
{
       int ii, jj, i, j, nkxz;
       int subx, subz;

       nkxz=nkx*nkz;

       fftw_complex* xin=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nkxz));
       fftw_complex* xout=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nkxz));

       fftw_plan xpi=fftw_plan_dft_2d(nkx,nkz,xin,xout,FFTW_BACKWARD,FFTW_ESTIMATE);

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

		xin[ii*nkz+jj]=xi[subx][subz]*I;
                //sf_warning("xi=%f xin=(%f, %f)",xi[subx][subz],creal(xin[ii*nkz+jj]),cimag(xin[ii*nkz+jj]));
	}

	fftw_execute(xpi);

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

		xo[subx][subz]=xout[ii*nkz+jj]/(nkxz);
	}

        fftw_destroy_plan(xpi);
        fftw_free(xin);
        fftw_free(xout);
}
