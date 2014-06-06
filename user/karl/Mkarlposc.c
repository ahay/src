/* n-D POCS interpolation using a general Lp-norm optimization
Note: Acquistion geometry represented by mask operator.
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <complex.h>
#include <fftw3.h>

#include "pthresh.h"

int main(int argc, char* argv[])
{
    	bool verb;
    	int i, i1, i2, index, n1, n2, num, dim, n[SF_MAX_DIM], nw, iter, niter, nthr;
	int npad;
	int npadded[SF_MAX_DIM],n1padded, n2padded, indxfile,ii[SF_MAX_DIM];
	int indxpadded;
    	float thr, pclip;
    	float *dobs_t, *thresh, *mask;
    	char key[7];
    	fftwf_complex *mm, *dd, *dobs;
    	fftwf_plan fft1, ifft1, fft2, ifft2;/* execute plan for FFT and IFFT */
    	sf_file in, out, Fmask;/* mask and I/O files*/ 


    	sf_init(argc,argv);	/* Madagascar initialization */

    	/* setup I/O files */
    	in=sf_input("in");	/* read the data to be interpolated */
    	out=sf_output("out"); 	/* output the reconstructed data */
    	Fmask=sf_input("mask");	/* read the (n-1)-D mask for n-D data */
 
    	if(!sf_getbool("verb",&verb))    	verb=false;
    	/* verbosity */
    	if (!sf_getint("niter",&niter)) 	niter=100;
    	/* total number iterations */
    	if (!sf_getfloat("pclip",&pclip)) 	pclip=10.;
    	/* starting data clip percentile (default is 99)*/
    	if (pclip <=0. || pclip > 100.)
	sf_error("pclip=%g should be > 0 and <= 100",pclip);

    	/* dimensions */
   	for (i=0; i < SF_MAX_DIM; i++) {
		snprintf(key,3,"n%d",i+1);
		if (!sf_getint(key,n+i) && 
		    (NULL == in || !sf_histint(in,key,n+i))) break;
		/*( n# size of #-th axis )*/  
		sf_putint(out,key,n[i]);
    	}

    	if (0==i) sf_error("Need n1=");
    	dim=i;

    	n1=n[0];
    	n2=sf_leftsize(in,1);
	for (i=0; i<SF_MAX_DIM; i++){
	  npadded[i]=1;
	}
	npad=10;
	/* if necessary increase n[0] to an even number, then pad */
	npadded[0]=((n[0]+1)/2)*2+npad; 
	n1padded=npadded[0];
	n2padded=1;
	for (i=1; i<dim; i++){
	  fprintf(stderr,"n[%d]=%d n2padded=%d\n",i,n[i],n2padded);
	  if(n[i]==1){
	    npadded[i]=n[i];
	    fprintf(stderr,"dim=%d and n[%d]=%d.  Why?\n",dim,i,n[i]);
	  }  else {
	    npadded[i]=n[i]+npad;
	  }
	  n2padded*=npadded[i];
	}
	
	nw=n1padded/2+1;
	num=nw*n2padded;/* total number of elements in frequency domain */
	fprintf(stderr,"n1padded=%d, n2padded=%d, nw=%d, num=%d\n",
		n1padded   , n2padded   , nw   , num);
    	/* allocate data and mask arrays */
	thresh=(float*)            malloc(nw      *n2padded*sizeof(float));
    	dobs_t=(float*)      fftwf_malloc(n1padded*n2padded*sizeof(float)); 
    	dobs=(fftwf_complex*)fftwf_malloc(nw      *n2padded*sizeof(fftwf_complex));
    	dd=(fftwf_complex*)  fftwf_malloc(nw      *n2padded*sizeof(fftwf_complex));
    	mm=(fftwf_complex*)  fftwf_malloc(nw      *n2padded*sizeof(fftwf_complex));

    	if (NULL != sf_getstring("mask")){
		mask=sf_floatalloc(n2padded);
    	} else sf_error("mask needed. the name of the mask file.");

	fprintf(stderr,"initialize data and mask\n");
	/* initialize the input data and mask arrays */
	memset(dobs_t,0,n1padded*n2padded*sizeof(float));
	memset(mask,0,n2padded*sizeof(float));
	for (indxfile=0; indxfile<n1*n2; indxfile+=n1){
	  sf_line2cart(dim,n,indxfile,ii);
	  indxpadded=sf_cart2line(dim,npadded,ii);
	  if(0)fprintf(stderr,
		  "indxpadded=%d,n1padded*n2padded=%d,indxpadded+n1=%d\n"
		  ,indxpadded   ,n1padded*n2padded   ,indxpadded+n1);
	  sf_floatread(&(dobs_t[indxpadded         ]),n1,in   );
	  sf_floatread(&(mask  [indxpadded/n1padded]), 1,Fmask);
	}

	fprintf(stderr,"fft1=fftwf_plan\n");
      /*fft1=fftwf_plan_many_dft_r2c(1, &n1      , n2      , dobs_t, &n1      , 1, n1      , dobs, &n1,       1, nw, FFTW_MEASURE);*/
     	fft1=fftwf_plan_many_dft_r2c(1, &n1padded, n2padded, dobs_t, &n1padded, 1, n1padded, dobs, &n1padded, 1, nw, FFTW_MEASURE);	
	fprintf(stderr,"ifft1=fftwf_plan\n");
   	ifft1=fftwf_plan_many_dft_c2r(1, &n1padded, n2padded, dobs, &n1padded, 1, nw, dobs_t, &n1padded, 1, n1padded, FFTW_MEASURE);
	fprintf(stderr,"fft2=fftwf_plan\n");
	fft2=fftwf_plan_many_dft(dim-1, &npadded[1], nw, mm, &npadded[1], nw, 1, mm, &npadded[1], nw, 1, FFTW_FORWARD,FFTW_MEASURE);
	fprintf(stderr,"ifft2=fftwf_plan\n");
	ifft2=fftwf_plan_many_dft(dim-1, &npadded[1], nw, mm, &npadded[1], nw, 1, mm, &npadded[1], nw, 1, FFTW_BACKWARD,FFTW_MEASURE);

	fprintf(stderr,"fftwf_execute\n");
	/*transform the data from time domain to frequency domain: tdat-->wdat*/
	fftwf_execute(fft1);
	fprintf(stderr,"fft normalize\n");
	for(i=0; i<num; i++) dobs[i]/=sqrtf(n1padded);
	fprintf(stderr,"initialize dd to 0\n");
	memset(dd,0,num*sizeof(fftwf_complex));

    	for(iter=0; iter<niter; iter++)
    	{
	        /* mm<--A^t dd */
	        fprintf(stderr,"memcpy\n");
		memcpy(mm, dd, num*sizeof(fftwf_complex));
		fprintf(stderr,"fftw\n");
		fftwf_execute(fft2);
		for(i=0; i<num; i++) mm[i]/=sqrtf(n2padded);

		fprintf(stderr,"threshold\n");
		/* perform hard thresholding: mm<--T{mm} */
		for(i=0; i<num; i++)	thresh[i]=cabsf(mm[i]);

	   	nthr = 0.5+num*(1.-0.01*pclip);
		fprintf(stderr,"num=%d,nthr=%d\n",num,nthr); 
	    	if (nthr < 0) nthr=0;
	    	if (nthr >= num) nthr=num-1;
		thr=sf_quantile(nthr,num,thresh);

		fprintf(stderr,"initial thr=%e\n",thr);
		thr*=((float)(niter-iter))/niter;
		fprintf(stderr,"thr=%e\n",thr);
		/* thr*=powf(0.01,(iter-1.0)/(niter-1.0)); */

		for(i=0; i<num; i++) mm[i]*=(cabsf(mm[i])>thr?1.:0.);

		/* mm<--A mm*/
		fftwf_execute(ifft2);

		for(i=0; i<num; i++) mm[i]/=sqrtf(n2padded);		

		/* d_rec = d_obs+(1-M)*A T{ At(d_rec) } */
		for(i2=0; i2<n2padded; i2++)
		for(i1=0; i1<nw; i1++)
		{ 
			index=i1+nw*i2;
			dd[index]=dobs[index]+(1.-mask[i2])*mm[index];
		}

		if (verb)    sf_warning("iteration %d;",iter+1);
    	}

	/*transform the data from frequency domain to time domain: wdat-->tdat*/
	memcpy(dobs, dd, num*sizeof(fftwf_complex));
	fftwf_execute(ifft1);
	for(i=0; i<n1padded*n2padded; i++) dobs_t[i]/=sqrtf(n1padded);

	/* write the data from the padded memory array to output file */	
    	for (indxfile=0; indxfile<n1*n2; indxfile+=n1){
	  sf_line2cart(dim,n,indxfile,ii);
	  indxpadded=sf_cart2line(dim,npadded,ii);
	  sf_floatwrite(&(dobs_t[indxpadded]),n1,out   );
	}

	free(thresh);
	fftwf_free(dobs_t);	
	fftwf_free(dobs);
	fftwf_free(mm);
	fftwf_free(dd);

    	exit(0);
}
