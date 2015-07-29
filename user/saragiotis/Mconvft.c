/* Trace-by-trace or data-by-trace convolution using Fourier transform. */
/* other can be a dataset with the same dimensions as in (except for axis) or a single trace.*/
/*
  Copyright (C) 2012 King Abdullah University of Science and Technology,
  Thuwal, Saudi Arabia
       
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

int main(int argc, char* argv[])
{
    int  i1,i2,i3,iw, n1,n2,n3,n12, m2,m12, nm2; 
    int  n[SF_MAX_DIM],m[SF_MAX_DIM], dim,mdim, axis, nfft,nw;
    float d2,dy2, o2,oy2, scale;
    float *x,*y=NULL,*z, *x1,*y1,*z1;
    char key1[3],key2[3];
    sf_complex *fft1,*fft2;
    kiss_fftr_cfg cfg,icfg;
    bool one=false;

    sf_file in, conv, other; 
    
    /* initialization */
    sf_init(argc,argv);
    in    = sf_input("in");
    other = sf_input("other");
    conv  = sf_output("out");
    
    
    /* Get rsf_file parameters */
    if (!sf_getint("axis",&axis)) axis=1;
    /*across which axis to convolve.*/
	
    dim = sf_filedims(in,n);
    if (dim<axis) sf_error("dim<axis (%d<%d).",dim,axis);
    n1 = n3 = 1;
    n2 = n[axis-1];
    for (i1=1; i1<=dim; i1++) {
	if (i1<axis)
	    n1 *=n[i1-1];
	else if(i1>axis)
	    n3 *= n[i1-1];
    }
	
    mdim = sf_filedims(other,m);
	
    if (mdim==1) {
	one = true;
	m2 = m[0];
    } else if (mdim!=dim) {
	sf_error("dimensions mismatch: dim of 'in' is %d but dim of 'other' is %d.",
		 dim,mdim);
	m2 = 0;
    } else {
	m2 = m[axis-1];
	for (i1=1; i1<=dim; i1++) {
	    if (i1!=axis && n[i1-1]!=m[i1-1])
		sf_error("dimension mismatch for axis %d.",i1);
	}
    }

    /* Check sampling consistency */
    snprintf(key1,3,"d%d",axis);
    snprintf(key2,3,"o%d",axis);
    if (!sf_histfloat(in, key1,&d2)) d2=1.;	
    if (!sf_histfloat(in, key2,&o2)) o2=0.;
	
    if (one) {
	if (!sf_histfloat(other,"d1",&dy2) || dy2!=d2) 
	    sf_error("%s of 'in' does not agree with d1 of 'other'.",key1);
	if (!sf_histfloat(other,"o1",&oy2)) oy2=0.;
    } else {
	if (!sf_histfloat(other,key1,&dy2) || dy2!=d2) 
	    sf_error("%s of 'in' and 'other' do not agree.",key1);
	if (!sf_histfloat(other,key2,&oy2)) oy2=0.;
    }
	
	
    nm2 = n2 + m2 - 1;
    snprintf(key1,3,"n%d",axis);	
    sf_putint  (conv,key1,nm2);
    sf_putfloat(conv,key2,o2+oy2);
	
    /* set up FFT */
    nfft  = 2*kiss_fft_next_fast_size((nm2+1)/2);
    nw    = nfft/2+1;
    /* dw    = 1./(nfft*d2); */
    scale = 1./nfft;
    cfg   = kiss_fftr_alloc(nfft,0,NULL,NULL);
    icfg  = kiss_fftr_alloc(nfft,1,NULL,NULL);
    fft1  = sf_complexalloc(nw);	
    fft2  = sf_complexalloc(nw);
    z  = sf_floatalloc(nfft*n1);
    z1 = sf_floatalloc(n1);
	
    n12 = n1*n2;
    m12 = one? m2:n1*m2;
    x  = sf_floatalloc(n12);
    if (!one) y = sf_floatalloc(m12);
    x1 = sf_floatalloc(nfft); for (i2=n2; i2<nfft; i2++) x1[i2] = 0.;
    y1 = sf_floatalloc(nfft); for (i2=m2; i2<nfft; i2++) y1[i2] = 0.;
	
    if (one) {
	sf_floatread(y1,m2,other);
	kiss_fftr(cfg, y1, (kiss_fft_cpx*) fft2);
    }
	
    for (i3=0; i3<n3; i3++) {
	sf_floatread(x, n12, in);
	if (!one) sf_floatread(y,m12,other);

	for (i1=0; i1<n1; i1++) {
	    for (i2=0; i2<n2; i2++)		/* transpose x */
		x1[i2] = x[i1+i2*n1];
	    kiss_fftr(cfg, x1, (kiss_fft_cpx*) fft1);

	    if (!one) {
		for (i2=0; i2<m2; i2++) /* transpose y */
		    y1[i2] = y[i1+i2*n1];
		kiss_fftr(cfg, y1, (kiss_fft_cpx*) fft2);
	    }
			
	    for (iw=0; iw<nw; iw++) {
#ifdef SF_HAS_COMPLEX_H
		fft1[iw] = fft1[iw]*fft2[iw]*scale;
#else
		fft1[iw] = sf_crmul(sf_cmul(fft1[iw],fft2[iw]), scale);
#endif
	    }
	    kiss_fftri(icfg,(kiss_fft_cpx*) fft1, z+i1*nfft);
	}
		
	for (iw=0; iw<nm2; iw++) {	/* transpose back */
	    for (i1=0; i1<n1; i1++)
		z1[i1] = z[iw+nfft*i1];
	    sf_floatwrite(z1,n1,conv);			
	}
	
    }


    exit(0);
}


