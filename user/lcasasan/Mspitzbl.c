/* Missing data interpolation in 2-D using F-X Prediction Error Filters
   based on Seismic trace interpolation in the F-X domain
   S.Spitz Geophysics 56, 785(1991). 
   This implementation is for bandlimited [f1 f2] signals

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
#include <math.h>
#include <float.h>
#include "estimLPF.h"
#include "predictionP.h"
#include "filtMCHF.h"


 #define b2(x)   (   (x) | (   (x) >> 1) )
 #define b4(x)   ( b2(x) | ( b2(x) >> 2) )
 #define b8(x)   ( b4(x) | ( b4(x) >> 4) )
 #define b16(x)  ( b8(x) | ( b8(x) >> 8) )  
 #define b32(x)  (b16(x) | (b16(x) >>16) )
 #define next_power_of_2(x)      (b32(x-1) + 1)	

void moment(float *data, int n, float *mean, float* var);

int main(int argc, char* argv[])
{
    float EPS=100*FLT_EPSILON;
	bool verb, inv, norm;    
    int i1,i2, n1, n2;    
    int order, ntraces;  /* input parameters */
    float o1, d1, o2, d2; /* data parameters */   
	float f[2];     /* bandwitch real and normalized f = fs * k */ 
    int   k[2];
	/* fft */
    int ny,nt,nyL,ntL,Ntot;	 
    kiss_fftr_cfg cfg,cfgL;	
	
	float *p,*pL,*p_out;
	float mean=0.0,var=0.0; /*for output normalization */
    kiss_fft_cpx *temp_fftL, *temp_fft, *out_fft;
	sf_complex *fftL, **FFTL, *fft,**FFT,**outFFT, **LPF;

    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT !=sf_gettype(in)) sf_error("Need float type");    

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
    
    //sf_warning("sampling frequency fs=%f",fs);

    if (!sf_getint("order",&order)) order=3;
    /* linear PEF order*/
    if (!sf_getint("ntraces",&ntraces)) ntraces=1;
    /* number of traces to be interpolated */
 	
	if (!sf_getfloat("f1",&f[0])) f[0] = 0.0;   /*lower  frequency in band limited signal >= 0.0  */
    if (!sf_getfloat("f2",&f[1])) f[1] = 1.0;   /*higher frequency in band limited signal <= 1.0  (normalized nyquist)*/
		//sf_warning("sampling frequency f[0]=%f f[1]=%f",f[0],f[1]);
   
	if ( (f[0] < 0.0) | (f[0] > (1.0) ) | (f[1] < 0.0) | (f[1] > (1.0)) ) 
		sf_error("### ERROR ### f[0]=%f and f[1]=%f must be between 0.0 and 1.0",f[0],f[1]);
	if (f[0]>= f[1])  
		sf_error("### ERROR ### f[0]=%f must be less f[1]=%f",f[0],f[1]);


	Ntot = (ntraces+1)*(n2-1)+1; /* total number of traces after interpolation */

 
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

	 if (!sf_getbool("norm",&norm)) norm = true;
    /* normalization flag */

    /* determine wavenumber sampling (for real to complex FFT) */
    //nt = 2*kiss_fft_next_fast_size((n1+1)/2);
	nt = next_power_of_2(n1);
    ntL = nt*(ntraces+1);

    ny = nt/2+1;
    nyL = ny*(ntraces+1);

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

	
    sf_putint  (out,"n2",Ntot);
	sf_putfloat  (out,"d2",d2/(ntraces+1));	

    for (i2=0; i2 < n2; i2++) {

	    sf_floatread (p,n1,in);
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
    
    /* normalized frequency [samples]*/
	k[0] = (int)(round ( f[0] * ny ) ) ;
	k[1] = (int)(round ( f[1] * ny ) ) ;	
	if(verb) sf_warning("sampled frequency ny=%d k[0]=%d k[1]=%d",ny,k[0],k[1]);
	
	if(verb) sf_warning("LPF start");
	/* estimation of linear prediction filter taps at each frequency in the input signal bandwitch*/
	LPFestim_init(n2,order,ntraces,false);
	for (i1=k[0]; i1<=k[1];i1++){
		LPFestim_apply(FFTL[i1], LPF[i1]); /* call LPF estimation */ 
	}
	LPFestim_close();
	if(verb) sf_warning("LPF end");
	
	/* we apply appendix B rules in Spitz paper to find 
	out of band LP filter components
	First we must create Q prediction filters which have components equal to 
	binomial coefficient(L,j) */

	if(verb) sf_warning("P prediction start");
	Pprediction_init(order, k, ny, false);

	Pprediction_apply(LPF);

	Pprediction_close();
	if(verb) sf_warning("P prediction end");

	/* application of multchannel filter to estimate missing traces */

	if(verb) sf_warning("MCHF start");
	MCHFfilt_init(n2,order,ntraces,false);

	for (i1=0; i1<ny-1;i1++){
		MCHFfilt_apply(LPF[i1], FFT[i1],outFFT[i1]); /* call MultiChannel Filtering */
	}

	MCHFfilt_close();
	if(verb) sf_warning("MCHF end");
	/* Inverse transform: we need to scale for 1/nt */;
	
	for (i2=0; i2 < Ntot; i2++){
		for (i1=0; i1<ny-1;i1++) {
			out_fft[i1].r = (1/(float) nt) * crealf( outFFT[i1][i2]);
			out_fft[i1].i = (1/(float) nt) * cimagf( outFFT[i1][i2]);
		}
		inv=true;
    	cfg =  kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
		kiss_fftri(cfg,out_fft,p_out);	


		if (norm) {
			moment(p_out, n1, &mean, &var);
			for (i1=0; i1 < n1; i1++)
//					p_out[i1]=( (p_out[i1]-mean)*sqrt(var) ) / ( sqrt(var) * ( sqrt(var) + EPS ) );
					p_out[i1]= (p_out[i1]-mean) / ( sqrt(var) + EPS );
		}
		sf_floatwrite((float*) p_out,n1,out);
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
	free(LPF[0]);     free(LPF);

    exit(0);
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

