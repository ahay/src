/* n-D IST interpolation using a generalized shrinkage operator and zero-padding
   Note: Acquistion geometry specified by mask operator
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)
  adapted from Mkarlpocs.c (Karl Schleicher)

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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "pthresh.h"

#ifdef SF_HAS_FFTW
#include <fftw3.h>

unsigned int nextpower2(unsigned int n)
/*< next power of 2 (>=n) >*/
{
    /* if n is a power of 2, return the same.  */
    if (!(n & (n-1)))     return (n);
 
    while (n & (n-1))    n = n & (n-1); 
    n = n << 1;
    return n;
}

int main(int argc, char* argv[])
{
    bool verb, pow2;
    char key[7], *mode;;
    int n1, n2, n1padded, n2padded, num, dim, n[SF_MAX_DIM], npadded[SF_MAX_DIM], ii[SF_MAX_DIM];
    int i, j, i1, i2, index, nw, iter, niter, nthr, *pad;
    float thr, pclip, normp;
    float *dobs_t, *thresh, *mask;
    fftwf_complex *mm, *dd, *dobs;
    fftwf_plan fft1, ifft1, fftrem, ifftrem;/* execute plan for FFT and IFFT */
    sf_file in, out, Fmask;	/* mask and I/O files*/ 

    sf_init(argc,argv);	/* Madagascar initialization */
    in=sf_input("in");	/* read the data to be interpolated */
    out=sf_output("out"); 	/* output the reconstructed data */
    Fmask=sf_input("mask");	/* read the (n-1)-D mask for n-D data */
 
    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity */
    if(!sf_getbool("pow2",&pow2))    	pow2=false;
    /* round up the length of each axis to be power of 2 */
    if (!sf_getint("niter",&niter)) 	niter=100;
    /* total number of iterations */
    if (!sf_getfloat("pclip",&pclip)) 	pclip=10.;
    /* starting data clip percentile (default is 10)*/
    if ( !(mode=sf_getstring("mode")) ) 	mode = "exp";
    /* thresholding mode: 'hard', 'soft','pthresh','exp';
       'hard', hard thresholding;	   'soft', soft thresholding; 
       'pthresh', generalized quasi-p;  'exp', exponential shrinkage */
    if (pclip <=0. || pclip > 100.)	sf_error("pclip=%g should be > 0 and <= 100",pclip);
    if (!sf_getfloat("normp",&normp)) 	normp=1.;
    /* quasi-norm: normp<2 */
    for (i=0; i < SF_MAX_DIM; i++) {/* dimensions */
	snprintf(key,3,"n%d",i+1);
	if (!sf_getint(key,n+i) && 
	    (NULL == in || !sf_histint(in,key,n+i))) break;
	/*( n# size of #-th axis )*/  
	sf_putint(out,key,n[i]);
    }
    if (0==i) sf_error("Need n1=");
    dim=i;
    pad=sf_intalloc (dim);
    for (i=0; i<dim; i++) pad[i]=0;
    sf_getints("pad",pad,dim); /* number of zeros to be padded for each axis */

    n1=n[0];
    n2=sf_leftsize(in,1);
    for (i=0; i<SF_MAX_DIM; i++) npadded[i]=1;
    npadded[0]=n1+pad[0];
    n1padded=npadded[0];
    n2padded=1;
    for (i=1; i<dim; i++){
	npadded[i]=n[i]+pad[i];
	if (pow2) {/* zero-padding to be power of 2 */
	    npadded[i]=nextpower2(n[i]);
	    fprintf(stderr,"n%d=%d n%dpadded=%d\n",i,n[i],i,npadded[i]);
	}
	n2padded*=npadded[i];
    }
    nw=npadded[0]/2+1;
    num=nw*n2padded;/* data: total number of elements in frequency domain */

    /* allocate data and mask arrays */
    thresh=(float*)            malloc(nw*n2padded*sizeof(float));
    dobs_t=(float*)      fftwf_malloc(n1padded*n2padded*sizeof(float));  /* time domain observation */
    dobs=(fftwf_complex*)fftwf_malloc(nw*n2padded*sizeof(fftwf_complex));/* freq-domain observation */
    dd=(fftwf_complex*)  fftwf_malloc(nw*n2padded*sizeof(fftwf_complex));
    mm=(fftwf_complex*)  fftwf_malloc(nw*n2padded*sizeof(fftwf_complex));
 
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2padded);
    } else sf_error("mask needed!");

    /* initialize the input data and mask arrays */
    memset(dobs_t,0,n1padded*n2padded*sizeof(float));
    memset(mask,0,n2padded*sizeof(float));
    for (i=0; i<n1*n2; i+=n1){
	sf_line2cart(dim,n,i,ii);
	j=sf_cart2line(dim,npadded,ii);
	sf_floatread(&dobs_t[j], n1, in);
	sf_floatread(&mask[j/n1padded], 1, Fmask);
    }

    /* FFT for the 1st dimension and the remaining dimensions */
    fft1=fftwf_plan_many_dft_r2c(1, &n1padded, n2padded, dobs_t, &n1padded, 1, n1padded, dobs, &n1padded, 1, nw, FFTW_MEASURE);
    ifft1=fftwf_plan_many_dft_c2r(1, &n1padded, n2padded, dobs, &n1padded, 1, nw, dobs_t, &n1padded, 1, n1padded, FFTW_MEASURE);
    fftrem=fftwf_plan_many_dft(dim-1, &npadded[1], nw, dd, &npadded[1], nw, 1, dd, &npadded[1], nw, 1, FFTW_FORWARD, FFTW_MEASURE);
    ifftrem=fftwf_plan_many_dft(dim-1, &npadded[1], nw, dd, &npadded[1], nw, 1, dd, &npadded[1], nw, 1, FFTW_BACKWARD, FFTW_MEASURE);

    /* transform the data from time domain to frequency domain: dobs_t-->dobs */
    fftwf_execute(fft1);
    for(i=0; i<num; i++) dobs[i]/=sqrtf(n1padded);
    memset(mm,0,num*sizeof(fftwf_complex));

    /* Iterative Shrinkage-Thresholding (IST) Algorithm:
       mm^{k+1}=T[mm^k+A^* M^* (dobs-M A mm^k)] (M^*=M; Mdobs=dobs)
       =T[mm^k+A^*(dobs-M A mm^k)]; (k=0,1,...niter-1)
       dd^=A mm^; 
    */
    for(iter=0; iter<niter; iter++)
    {
	/* dd<-- A mm^k */
	memcpy(dd, mm, num*sizeof(fftwf_complex));
	fftwf_execute(ifftrem);
	for(i=0; i<num; i++) dd[i]/=sqrtf(n2padded);

	/* apply mask: dd<--dobs-M A mm^k=dobs-M dd */
	for(i2=0; i2<n2padded; i2++)
	    for(i1=0; i1<nw; i1++)
	    { 
		index=i1+nw*i2;
		dd[index]=dobs[index]-mask[i2]*dd[index];
	    }

	/* mm^k += A^*(dobs-M A mm^k); dd=dobs-M A mm^k */
	fftwf_execute(fftrem);
	for(i=0; i<num; i++) mm[i]+=dd[i]/sqrtf(n2padded);		

	/* perform thresholding */
	for(i=0; i<num; i++)	thresh[i]=cabsf(mm[i]);
	nthr = 0.5+num*(1.-0.01*pclip);  /* round off */
	if (nthr < 0) nthr=0;
	if (nthr >= num) nthr=num-1;
	thr=sf_quantile(nthr, num, thresh);
	sf_cpthresh(mm, num, thr, normp, mode);

	if (verb) sf_warning("iteration %d;",iter+1);
    }

    /* frequency--> time domain: dobs-->dobs_t */
    memcpy(dd, mm, num*sizeof(fftwf_complex));
    fftwf_execute(ifftrem);
    for(i=0; i<num; i++) dd[i]/=sqrtf(n2padded);
    memcpy(dobs, dd, num*sizeof(fftwf_complex));
    fftwf_execute(ifft1);
    for(i=0; i<n1padded*n2padded; i++) dobs_t[i]/=sqrtf(n1padded);
	
    for (i=0; i<n1*n2; i+=n1){
	sf_line2cart(dim,n,i,ii);
	j=sf_cart2line(dim,npadded,ii);
	sf_floatwrite(&dobs_t[j],n1,out);
    }

    free(thresh);
    fftwf_free(dobs_t);
    fftwf_free(dobs);
    fftwf_free(dd);
    fftwf_free(mm);

    exit(0);
}
#endif
