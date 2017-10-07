/* Missing data interpolation in 2-D using F-X Prediction Error Filters

   Based on Seismic trace interpolation in the F-X domain
   S.Spitz Geophysics 56, 785(1991). 
*/
/*
  Copyright (C) 2010 Politecnico di Milano
  
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
#include <float.h>
#include "estimLPF.h"
#include "predictionP.h"
#include "filtMCHF.h"

#include "spitzbl.h"

#define b2(x)   (   (x) | (   (x) >> 1) )
#define b4(x)   ( b2(x) | ( b2(x) >> 2) )
#define b8(x)   ( b4(x) | ( b4(x) >> 4) )
#define b16(x)  ( b8(x) | ( b8(x) >> 8) )  
#define b32(x)  (b16(x) | (b16(x) >>16) )
#define next_power_of_2(x)      (b32(x-1) + 1)


void moment(float *data, int n, float *mean, float* var);


void spitzbl(float *in,float *out, int *n, int order, 
	     int ntraces, float* f, bool norm)
/*< main action >*/
{
    float EPS=100*FLT_EPSILON;
    bool verb=false, inv;    
    int i1,i2, n1, n2;       
    /* fft */
    int ny,nt,nyL,ntL,Ntot;
    int k[2]; /* frequency bandwitch normalized in samples*/
    float  *p,*pL,*p_out;
    kiss_fft_cpx *temp_fftL, *temp_fft, *out_fft;
    sf_complex *fftL, **FFTL, *fft,**FFT,**outFFT, **LPF;
    kiss_fftr_cfg cfg,cfgL;	
    float mean=0.0,var=0.0;	

    n1 = n[0];
    n2 = n[1];

    /* n12 = n1*n2; */
    //sf_warning("SONO DENTRO SPITZ.C n12 = %d",n12);
    
 	
    Ntot = (ntraces+1)*(n2-1)+1; /* total number of traces after interpolation */
    /* Nint = Ntot-n2; 			number of new traces */
    //sf_warning("total number of traces %d",Ntot);
    //sf_warning("total number of interp traces %d",Nint);
 
   
	
    /* determine wavenumber sampling (for real to complex FFT) */
    //nt = 2*kiss_fft_next_fast_size((n1+1)/2);
    nt = next_power_of_2(n1);
    /*sf_warning("\n fft samples = %d",nt);*/

    ntL = nt*(ntraces+1);
    /*sf_warning("\n fftL samples = %d",ntL);*/
    ny = nt/2+1;
    /*sf_warning("\n Nyquist sample = %d",ny);*/

    nyL = ny*(ntraces+1);
    /*sf_warning("\n Nyquist nyL sample = %d",nyL);*/
    

    p = sf_floatalloc(nt);
    pL = sf_floatalloc(ntL);
    p_out = sf_floatalloc(nt);


    temp_fftL = (kiss_fft_cpx*) sf_complexalloc(nyL);
    temp_fft  = (kiss_fft_cpx*) sf_complexalloc(ny);
    out_fft  = (kiss_fft_cpx*) sf_complexalloc(ny);

    for (i1=0; i1<ny; i1++) {
	out_fft[i1].r = 0.0;
	out_fft[i1].i = 0.0;
    }


    fftL 	= sf_complexalloc(nyL);
    fft  	= sf_complexalloc(ny);

    FFTL 	= sf_complexalloc2(n2,nyL);
    FFT 	= sf_complexalloc2(n2,ny);
    outFFT 	= sf_complexalloc2(Ntot,ny);

    LPF		= sf_complexalloc2(order,ny);
    
    inv=false;
    cfg =  kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
    cfgL = kiss_fftr_alloc(ntL,inv?1:0,NULL,NULL);    

	


    for (i2=0; i2 < n2; i2++) {

	//sf_floatread (p,n1,in);
	//sf_floatread (p,n1,in);
	for (i1=0;i1 < n1;i1++) {
	    p[i1] = in[i1+i2*n1]; 
	}
	for (i1=0; i1 < n1; i1++)
	    pL[i1]=p[i1];		

	for (i1=n1; i1 < nt; i1++) {  // Zero padding up to nt 
	    p[i1]=0.0;
	}
	    
        for (i1=n1; i1 < ntL; i1++) {  // Zero padding up to ntL 
	    pL[i1]=0.0;
	}
	    
	kiss_fftr (cfg,p,temp_fft);
	kiss_fftr (cfgL,pL,temp_fftL);
	    
	// cast to sf_complex to be consistent with the RSF API!
	    
        for (i1=0; i1 < ny; i1++) {
	    FFT[i1][i2]
		=sf_cmplx(temp_fft[i1].r , temp_fft[i1].i);
	}

	for (i1=0; i1 < nyL; i1++) {
	    FFTL[i1][i2]
		=sf_cmplx(temp_fftL[i1].r , temp_fftL[i1].i);	
	}


    }
    
    k[0] = (int)(round ( f[0]  * ny)) ;
    k[1] = (int)(round ( f[1]  * ny)) ;	
    //sf_warning("sampled frequency ny=%d f[0]=%f f[1]=%f k[0]=%d k[1]=%d",ny,f[0],f[1],k[0],k[1]);
	
    LPFestim_init(n2,order,ntraces,verb);
    for (i1=k[0]; i1<=k[1];i1++){
	LPFestim_apply(FFTL[i1], LPF[i1]); /* call LPF estimation */ 
    }

    LPFestim_close();
	
	
    /* we apply appendix B rules in Spitz paper to find 
       out of band LP filter components
       First we must create Q prediction filters which have components equal to 
       binomial coefficient(L,j) */

    Pprediction_init(order, k, ny, verb);

    Pprediction_apply(LPF);

    Pprediction_close();
	

    /* application of multchanell filter to estimate missing traces */

    MCHFfilt_init(n2,order,ntraces,verb);

    for (i1=0; i1<ny-1;i1++){
	MCHFfilt_apply(LPF[i1], FFT[i1],outFFT[i1]); /* call MultiCHannel Filtering */ 
    }

    MCHFfilt_close();
		
    // Inverse transform: we need to scale for 1/nt;
	
    for (i2=0; i2 < Ntot; i2++){
	for (i1=0; i1<ny-1; i1++) {
	    out_fft[i1].r = (1/(float) nt) * crealf( outFFT[i1][i2]);
	    out_fft[i1].i = (1/(float) nt) * cimagf( outFFT[i1][i2]);
	}
	inv=true;
	cfg =  kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
	kiss_fftri(cfg,out_fft,p_out);	

	if (norm) moment(p_out, n1, &mean, &var);

	for (i1=0; i1 < n1; i1++) {
	    if (norm)
//				out[i1+i2*n1]=( (p_out[i1]-mean)*sqrt(var) ) / ( sqrt(var) * ( sqrt(var) + EPS ) );
		out[i1+i2*n1]= (p_out[i1]-mean) / ( sqrt(var) + EPS );
	    else
		out[i1+i2*n1]=p_out[i1];

	}
    }	

    
    
    free(p); 
    free(pL); 
    free(p_out); 

    free(temp_fftL); 
    free(temp_fft); 
    free(out_fft);  
        
    free(fftL);
    free(fft);

    free(FFTL[0]);	  free(FFTL);
    free(FFT[0]);  	  free(FFT);
    free(outFFT[0]);  free(outFFT);

    free(LPF[0]);  free(LPF);
	 	
}


void moment(float *data, int n, float *mean, float *var) {
/*compute mean and variance*/
    int i;
	
    *mean=0.0;
    *var=0.0;
    if (n<=1) sf_error("n must be at least 2 in moment function");
    for (i=0;i<n;i++) {
	*mean += data[i];
	*var  +=data[i]*data[i];
    }
    *mean/=(float)n;
    //var=	[ sum(data^2) - n*mean^2 ] / [ n-1 ]
    *var=( *var-( (float)n * (*mean) * (*mean) ) ) / ( (float)(n-1) ) ;

}
