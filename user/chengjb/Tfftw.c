/* Test FFT code from FFTW3 */
/*
  Copyright (C) 2012 Tongji University (Jiubing Cheng) 
  and The University of Texas at Austin (Sergey Fomel)
 
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
#include <fftw3.h>

/*****************************************************************************************/
int main(int argc, char* argv[])
{

       int   m=11, i;

       sf_warning("point number m=%d",m);
       sf_warning("========================================");


#ifdef SF_HAS_FFTW  /* using FFTW in Madagascar */
 
       sf_complex *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(m);

       xp=fftwf_plan_dft_1d(m, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_1d(m,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: forward */
       for(i=0;i<m;i++) 
       {
          xin[i]=sf_cmplx(sin(i*1.0), 0.);
          sf_warning("i=%d xin=(%f,%f)",i,creal(xin[i]));
       }
       fftwf_execute(xp);
           
       sf_warning("========================================");
       for(i=0;i<m;i++) 
         sf_warning("i=%d xout=(%f,%f)",i,creal(xout[i]),cimag(xout[i]));

       for(i=0;i<m;i++) 
          xin[i]=xout[i];

       fftwf_execute(xpi);
           
       sf_warning("========================================");
       for(i=0;i<m;i++) 
         sf_warning("i=%d xout=(%f,%f)",i,creal(xout[i])/m,cimag(xout[i])/m);

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
       free(xin);
       free(xout);
#else  /* using FFTW in user's own computer */


       fftw_complex *xin, *xout;

       fftw_plan xp;
       fftw_plan xpi;

       xin=fftw_complexalloc(m);
       xout=fftw_complexalloc(m);

       xp=fftw_plan_dft_1d(m, (fftw_complex *) xin, (fftw_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftw_plan_dft_1d(m,(fftw_complex *) xin, (fftw_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: forward */
       for(i=0;i<m;i++) 
          xin[i]=sf_cmplx(1.0, 0.);

       fftw_execute(xp);
           
       for(i=0;i<m;i++) 
          sf_warning("i=%d xout=(%f,%f)",i,creal(xout[i]),cimag(xout[i]));

       for(i=0;i<m;i++) 
          xin[i]=xout[i];

       fftw_execute(xpi);
           
       for(i=0;i<m;i++) 
         sf_warning("i=%d xout=(%f,%f)",i,creal(xout[i])/m,cimag(xout[i])/m);

       fftw_destroy_plan(xp);
       fftw_destroy_plan(xpi);
       free(xin);
       free(xout);
#endif
       exit(0);
}
