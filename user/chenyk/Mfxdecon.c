/* Random noise attenuation using f-x deconvolution */
/*
  Copyright (C) 2013 University of Texas at Austin
   
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

  References:      							
  Canales(1984):'Random noise reduction' 54th. SEGM	
  Gulunay(1986):'FXDECON and complex Wiener Predicition filter' 56th. SEGM	                
  Galbraith(1991):'Random noise attenuation by F-X prediction: a tutorial' 61th. SEGM
  Seismic Unix (SU) software package: http://www.cwp.mines.edu/cwpcodes/.
*/
#include <rsf.h>
#include "fxdeconutil.h"


int main(int argc, char *argv[])
{
    int n1;			/* number of samples in input trace	 */
    int n1w;		/* number of samples in time window      */
    int N1w;	       	/* number of time windows		*/
    int n1ws;	       	/* number of samples in window (wo taper)*/
    int n1wf;	       	/* number of samples in first window     */
    int n1wi;	       	/* number of samples in intermed. window */
    int n1taper;	       	/* number of samples in taper		*/
    int n2;	       		/* number of input traces		*/
    int n2w;	       	/* number of traces in window		*/
    int N2w;	       	/* number of spacial windows		*/
    int n2wu;	       	/* updated number of traces in window	*/
    int nfft;		/* transform length			*/
    int nf;			/* number of frequencies		*/
    int lenf;		/* number of traces for filter		*/
    int i2w,i2,jr,itt;	/* summation indexes			*/
    int i1,i1w,ifq,ir;	/* summation indexes			*/
    int ig,ifv;		/* summation indexes			*/
    bool verb;		/* flag to get advisory messages	*/
    int *ipvt;		/* indices of pivot permutations	*/

    float *info;		/* index of last zero pivot		*/
    float dt;		/* sampling interval in secs		*/
    float df;		/* sample interval in Hz		*/
    float fmin;		/* minimum frequency to process in Hz	*/
    float fmax;		/* maximum frequency to process in Hz	*/
    float taper;		/* length of taper			*/
    float twlen;       	/* time window length 			*/
    float **tracein;      	/* real trace			   	*/
    float **traceout;     	/* real trace			   	*/
    float *traceintw;       /* real trace in time window - input   	*/
    float *traceinttw;      /* real trace in time window - input   	*/
    float *traceotw;       	/* real trace in time window - output   */
    float *traceottw;      	/* real trace in time window - output   */
    float *rauto;       	/* real part of autocorrelation   	*/
    float *iauto;       	/* imaginary part of autocorretation 	*/
    float **cautom;      	/* complex autocorrelation matrix	*/
    float *aav;       	/* advanced autocorrelation vector	*/

    kiss_fft_cpx *cdata;	/* complex transformed trace (window)  	*/
    kiss_fft_cpx **fdata;	/* data - freq. domain          	*/
    kiss_fft_cpx **fdataw;	/* data - freq. domain - in window     	*/
    kiss_fft_cpx *sfv;	/* single frequency vector 	  	*/
    kiss_fft_cpx **ffv;	/* filtered frequency vector 	  	*/
    kiss_fft_cpx *fv;	/* frequency vector		      	*/
    kiss_fft_cpx *sfvout;	/* single frequency output vector      	*/
    kiss_fft_cpx *sfvcj;	/* complex conjugate of sfv	  	*/
    kiss_fft_cpx *acorr;	/* autocorrelation output	  	*/
    kiss_fft_cpx *filter;	/* filter               	  	*/

    kiss_fftr_cfg forw, invs;
    sf_init(argc,argv);
    sf_file in, out;

    in = sf_input("in");
    out = sf_output("out");

    if(!sf_getbool("verb",&verb)) verb=false;
    /* flag to get advisory messages */
	
    if(!sf_histint(in, "n1", &n1))  sf_error("No n1 in input");  if (verb) sf_warning("n1 = %i",n1);
    if(!sf_histfloat(in, "d1", &dt)) sf_error("No d1 in input"); if (verb) sf_warning("dt= %f",dt);
    if(!sf_histint(in,"n2",&n2))   sf_error("No n2 in input");   if (verb) sf_warning("n2= %f",n2);


    if(!sf_getfloat("taper",&taper)) taper=.1;
    /* length of taper */
    if (taper==0.0) taper=.004; 

    if(!sf_getfloat("fmin",&fmin)) fmin=1.;    
    /* minimum frequency to process in Hz */
    if (fmin==0.0)  if (verb) sf_warning("using fmin=1 Hz");

    if(!sf_getfloat("fmax",&fmax)) fmax=1./(2*dt); 
    /* maximum frequency to process in Hz */
    if (fmax==0.0) if (verb) sf_warning("using fmax=1/(2*dt) Hz");

    if (!sf_getfloat("twlen", &twlen)) twlen=(float)(n1-1)*dt;
    /* time window length */
    if (twlen<.3) {
	twlen=.3; if (verb) sf_warning("twlen cannot be less than .3s, using .3s");
    }
    /* setting taper and spatial and temporal windows */
    n1taper =round(taper/dt);
    n1taper = (n1taper%2 ? n1taper+1 : n1taper); 

    N1w = round((n1-1)*dt/twlen);

    if (N1w==1) taper=0.0;

    n1ws = round(twlen/dt) + 1;
    n1wf = n1ws + n1taper/2;
    n1wi = n1ws + n1taper;

    if (!sf_getint("n2w", &n2w)) n2w=10;
    /* number of traces in window */
    if (!n2w) {
	n2w = 10;
	if (verb) sf_warning("n2w cannot be zero, using 10 traces");
    }	
    if (verb) sf_warning("n2w = %i",n2w);

    n2wu = n2w;

    if (!sf_getint("lenf", &lenf)) lenf=4; 
    /* number of traces for filter */
    if (!lenf) if (verb) sf_warning("using lenf=4");	
	
    N2w = n2/n2w;

    /* Computute FFT optimization number */
    nfft = 2*kiss_fft_next_fast_size((n1wf+1)/2);

    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);

    nf = nfft/2 + 1;
    df = 1.0/(nfft*dt);

    /* space allocation */
    cdata 	= (kiss_fft_cpx*) sf_complexalloc(nfft);
    traceintw =sf_floatalloc(nfft);
    traceinttw=sf_floatalloc(nfft);
    fdata   =(kiss_fft_cpx**) sf_complexalloc2(nf,n2);
    fdataw  =(kiss_fft_cpx**) sf_complexalloc2(nf,2*n2w+2*lenf);
    tracein = sf_floatalloc2(n1,n2);
    traceout= sf_floatalloc2(n1,n2);
    sfv   = (kiss_fft_cpx*) sf_complexalloc(2*n2w+2*lenf);
    sfvcj = (kiss_fft_cpx*) sf_complexalloc(2*n2w+2*lenf);
    acorr= (kiss_fft_cpx*) sf_complexalloc(lenf+1);
    rauto  = sf_floatalloc(lenf+1);
    iauto  = sf_floatalloc(lenf+1);
    cautom = sf_floatalloc2(2*lenf,2*lenf);
    aav = sf_floatalloc(2*lenf);
    filter = (kiss_fft_cpx*) sf_complexalloc(2*lenf+1);
    ffv   = (kiss_fft_cpx**) sf_complexalloc2(nf,2*n2w);
    sfvout=(kiss_fft_cpx*) sf_complexalloc(2*n2w);
    fv   =(kiss_fft_cpx*) sf_complexalloc(nfft);
    traceotw = sf_floatalloc(nfft);
    traceottw= sf_floatalloc(nfft);
    ipvt 	= sf_intalloc(4*n2w);
    info	= sf_floatalloc(4*n2w);

    /* zero output file */
    memset((void *) traceout[0], 0, n2*n1*sizeof(float));

    /* load traces into the zero-offset array and close tmpfile */
    sf_floatread(tracein[0],n1*n2,in);
	
	
    /* If dt not set, issue advisory on frequency step df */
    if (dt && verb) sf_warning("df=%f", 1.0/(nfft*dt));

    if (verb) sf_warning("nf=%i, df=%f, nfft=%i, n1taper=%i", nf,df,nfft,n1taper);

    /* loop over time windows */
    for (i1w=0;i1w<N1w;i1w++) {

	if (i1w>0 && i1w<N1w-1) n1w=n1wi; 
	else if (i1w==0) 
	    if (N1w>1) n1w = n1wf;
	    else        n1w = n1;
	else
	    n1w = n1 - n1ws*i1w + n1taper/2;
 
	if (verb) sf_warning("i1w=%i, N1w=%i, n1w=%i, twlen=%f", i1w,N1w,n1w,twlen); 

	/* zero fdata */
	memset((void *) fdata[0], 0, nf*n2*sizeof(kiss_fft_cpx));
     
	/* select data */
	for (i2=0;i2<n2;i2++) {
 
	    if (i1w>0)
		for (i1=0;i1<n1w;i1++)
		    traceintw[i1]=tracein[i2][i1 + i1w*n1ws - n1taper/2];
	    else
		for (i1=0;i1<n1w;i1++)  
		    traceintw[i1]=tracein[i2][i1];	  

	    memset((void *) (traceintw + n1w), 0, (nfft-n1w)*sizeof(float));
	    memset((void *) cdata, 0, nfft*sizeof(kiss_fft_cpx));

	    /* FFT from t to f */
	    for (i1=0;i1<nfft;i1++)
		traceinttw[i1]=(i1%2 ? -traceintw[i1] : traceintw[i1]);
	    kiss_fftr(forw, traceinttw, cdata);

	    /* Store values */    
	    for (ifq = 0; ifq < nf; ifq++) { 
		fdata[i2][ifq] = cdata[nf-1-ifq];
	    }
	}

	/* Loop over space windows */
	for (i2w=0;i2w<N2w;i2w++){

	    /* to take care of a possible incomplete last window */
	    if (n2<i2w*n2w+2*n2w) 
		n2wu = n2 - i2w*n2w;
	    else
		n2wu = n2w;

	    if (verb) {
		sf_warning("i2w=%i, n2=%i, n2w=%i",
			   i2w,n2,n2w);
		sf_warning("n2wu=%i, N2w=%i, lenf=%i",
			   n2wu,N2w,lenf);
	    }

	    /* zero fdataw */
	    for (i2=0;i2<n2w+2*lenf;i2++)
		memset((void *) fdataw[i2], 0, nf*sizeof(kiss_fft_cpx));

	    /* select data */
	    for (i2=0;i2<n2wu+2*lenf;i2++) 
		for (ifq = 0; ifq < nf; ifq++) {

		    if (i2w>0 && i2w<N2w-1)  
			fdataw[i2][ifq] = fdata[i2 + i2w*n2w - lenf][ifq];
		    else if (i2w==0)
			if (i2>=lenf && i2<n2w+lenf) 
			    fdataw[i2][ifq] = fdata[i2 - lenf][ifq];
			else if (i2<lenf) 
			    fdataw[i2][ifq] = fdata[0][ifq];
			else 
			    if (N2w>1) 
				fdataw[i2][ifq] = fdata[i2 - lenf][ifq];
			    else 
				fdataw[i2][ifq] = fdata[n2-1][ifq];
		    else
			if (i2<n2wu+lenf)
			    fdataw[i2][ifq] = fdata[i2 + i2w*n2w - lenf][ifq];
			else 
			    fdataw[i2][ifq] = fdata[n2-1][ifq];
		}

	    /* loop over frequencies */
	    for (ifq=0;ifq<nf;ifq++) {

		if ((float)ifq*df>=fmin && (float)ifq*df<=fmax) {

		    /* Loop over space window */
		    memset((void *) sfv, 0, (n2wu+2*lenf)*sizeof(kiss_fft_cpx));
		    memset((void *) sfvcj, 0, (n2wu+2*lenf)*sizeof(kiss_fft_cpx));

		    for (i2=0;i2<n2wu+2*lenf;i2++) {
	  
			sfv[i2]=fdataw[i2][ifq];
			sfvcj[i2]=sf_conjf(fdataw[i2][ifq]);
		    }

		    memset((void *) acorr, 0, (lenf+1)*sizeof(kiss_fft_cpx));

		    /* complex autocorrelation */
		    cxcor(n2wu,0,sfv,n2wu,0,sfv,lenf+1,0,acorr);

		    /* zeroing files */
		    memset((void *) rauto, 0, (lenf+1)*sizeof(float));
		    memset((void *) iauto, 0, (lenf+1)*sizeof(float));

		    /* taking real and imaginary parts */
		    for (i2=0;i2<lenf+1;i2++) {
			rauto[i2]=acorr[i2].r;
			iauto[i2]=acorr[i2].i;
		    }

		    /* zeroing files */
		    memset((void *) aav, 0, 2*lenf*sizeof(float));
		    memset((void *) filter, 0, (2*lenf+1)*sizeof(kiss_fft_cpx));
		    for (ir=0;ir<2*lenf;ir++) 
			memset((void *) cautom[ir], 0, 2*lenf*sizeof(float));

		    /* matrix problem */
		    for (ir=0;ir<lenf;ir++) 
			for (jr=0;jr<lenf;jr++) { 
			    if (ir>=jr) cautom[ir][jr]=acorr[ir-jr].r;
			    else        cautom[ir][jr]=acorr[jr-ir].r;
			}

		    for (ir=lenf;ir<2*lenf;ir++)
			for (jr=0;jr<lenf;jr++) {
			    if (ir-lenf<jr) cautom[ir][jr]=-acorr[jr-ir+lenf].i;
			    else            cautom[ir][jr]= acorr[ir-jr-lenf].i;
			}

		    for (ir=lenf;ir<2*lenf;ir++)
			for (jr=lenf;jr<2*lenf;jr++)
			    cautom[ir][jr]=cautom[ir-lenf][jr-lenf];

		    for (ir=0;ir<lenf;ir++)
			for (jr=lenf;jr<2*lenf;jr++)
			    cautom[ir][jr]=-cautom[ir+lenf][jr-lenf];

		    for (ig=0;ig<2*lenf;ig++) {
			if (ig<lenf) aav[ig]=acorr[ig+1].r;
			else aav[ig]=acorr[ig-lenf+1].i;
		    }

		    lu_decomposition(2*lenf,cautom,ipvt,info);
		    backward_substitution(2*lenf,cautom,ipvt,aav);
      
		    /* construct filter */
		    for (ifv=0,ig=lenf-1;ifv<lenf;ifv++,ig--) 
			filter[ifv]=sf_conjf(cmplx(aav[ig]/2.,aav[ig+lenf]/2.));

		    for (ifv=lenf+1,ig=0;ifv<2*lenf+1;ifv++,ig++) 
			filter[ifv]=cmplx(aav[ig]/2.,aav[ig+lenf]/2.);
	 
		    memset((void *) sfvout, 0, n2wu*sizeof(kiss_fft_cpx));

		    /* convolution of data with filter */
		    /* output is one sample ahead */
		    cconv(n2wu+2*lenf,-lenf,sfv,2*lenf+1,-lenf,filter,n2wu,0,sfvout); 

		    /* store filtered values */
		    for (i2=0;i2<n2wu;i2++) ffv[i2][ifq]=sfvout[i2];

		}
	    } /* end of frequencies loop */

	    /* loop along space windows */
	    for (i2=0;i2<n2wu;i2++) {
    
		/* select data */
		for (ifq=0,itt=nf-1;ifq<nf;ifq++,itt--)
		    fv[ifq] = ffv[i2][itt]; 

		memset((void *) (fv+nf), 0, (nfft-nf)*sizeof(kiss_fft_cpx));
		memset((void *) traceotw, 0, nfft*sizeof(float));

		/* FFT back from f to t and scaling */
		kiss_fftri(invs,fv,traceotw);
		for (i1=0;i1<SF_MIN(n1,nfft);i1++)
		    traceotw[i1]/=nfft; 
		for (i1=0;i1<SF_MIN(n1,nfft);i1++)
		    traceottw[i1]=(i1%2 ? -traceotw[i1] : traceotw[i1]); 
      
		/*loop along time */
		if (N1w>1) {
		    /* first portion of time window */
		    if (i1w>0) 
			for (i1=0;i1<n1taper;i1++)
			    traceout[i2w*n2w+i2][i1+i1w*n1ws-n1taper/2]+=
				traceottw[i1]*((float)(i1)*dt/taper);
		    else 
			for (i1=0;i1<n1taper;i1++)
			    traceout[i2w*n2w+i2][i1]=traceottw[i1];

		    /* intermediate portion of time window */
		    if (i1w>0) 
			for (i1=n1taper;i1<n1w-n1taper;i1++)
			    traceout[i2w*n2w+i2][i1+i1w*n1ws-n1taper/2]=traceottw[i1];
		    else 
			for (i1=n1taper;i1<n1w-n1taper;i1++)
			    traceout[i2w*n2w+i2][i1]=traceottw[i1];

		    /* last portion of time window */
		    if (i1w>0 && i1w<N1w-1) 
			for (i1=n1w-n1taper;i1<n1w;i1++)
			    traceout[i2w*n2w+i2][i1+i1w*n1ws-n1taper/2]+=
				traceottw[i1]*(1.-((float)(i1-n1w+n1taper))*dt/taper);
		    else if (i1w==N1w-1)
			for (i1=n1w-n1taper;i1<n1w;i1++)
			    traceout[i2w*n2w+i2][i1+i1w*n1ws-n1taper/2]=traceottw[i1];
		    else 
			for (i1=n1w-n1taper;i1<n1w;i1++)
			    traceout[i2w*n2w+i2][i1]+=traceottw[i1]*(1.-((float)(i1-n1w+n1taper))*dt/taper);
		}
		else {
		    for (i1=0;i1<n1;i1++) 
			traceout[i2w*n2w+i2][i1]=traceottw[i1];
		}

	    } /* end loop over space windows */

	} /* end loop over space windows */

    } /* end of time windows loop */

    /* Free allocated memory */
    free(traceintw);
    free(traceinttw);
    free(fdataw[0]);
    free(fdataw);
    free(sfv);
    free(sfvcj);
    free(acorr);
    free(rauto);
    free(iauto);
    free(cautom[0]);
    free(cautom);
    free(aav);
    free(filter);
    free(sfvout);
    free(fv);
    free(traceotw);
    free(traceottw);
    free(ffv[0]);
    free(ffv);
    free(tracein[0]);
    free(tracein);
 
    /* Write output data to file */
    sf_floatwrite(traceout[0], n1*n2, out);

    exit(0);
}







