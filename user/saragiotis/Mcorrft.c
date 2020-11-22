/* Trace-by-trace or data-by-trace correlation using Fourier transform. 
other can be a dataset with the same dimensions as in (except for axis) or a single trace.
 If other is not specified, auto-correlation is assumed.*/
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
    int  i1,i2,i3,iw, n1,n2,n3,n12, m2=0,m12=0, nm2; 
	int  n[SF_MAX_DIM],m[SF_MAX_DIM], dim,mdim, axis, nfft,nw;
    float d2,dy2, o2,oy2, scale,shift;
    float *x,*x1, *y=NULL,*y1=NULL, *z,*z1, tmp;
	char key[3];
	sf_complex *fft1,*fft2=NULL, *ce=NULL;
	kiss_fftr_cfg cfg,icfg;
	bool one=false, autocorr=false;

    sf_file in, corr, other=NULL; 
    
    // initialization
    sf_init(argc,argv);
    in    = sf_input("in");
    corr  = sf_output("out");
    
    
    // Get rsf_file parameters */
    if (!sf_getint("axis",&axis)) axis=1;
	/*across which axis to correlate.*/
	snprintf(key,3,"d%d",axis);  if (!sf_histfloat(in, key, &d2)) d2=1.;	
	snprintf(key,3,"o%d",axis);  if (!sf_histfloat(in, key, &o2)) o2=0.;
	
	
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
	
	
	if (NULL!=(sf_getstring("other"))) {
		other = sf_input("other");
		mdim = sf_filedims(other,m);
	
		if (mdim==1) {
			one = true;
			m2 = m[0];
		} else if (mdim!=dim) {
			sf_error("dimensions mismatch: dim of 'in' is %d but dim of 'other' is %d.",
					 dim,mdim);
		} else {
			m2 = m[axis-1];
			for (i1=1; i1<=dim; i1++) {
				if (i1!=axis && n[i1-1]!=m[i1-1])
					sf_error("dimension mismatch for axis %d.",i1);
			}
		}
		
		// Check sampling consistency	
		snprintf(key,3,"d%d",axis);
		if (one) {
			if (!sf_histfloat(other,"d1",&dy2) || dy2!=d2) 
				sf_error("%s of 'in' does not agree with d1 of 'other'.",key);
			if (!sf_histfloat(other,"o1",&oy2)) oy2=0.;
		} else {
			if (!sf_histfloat(other,key,&dy2) || dy2!=d2) 
				sf_error("%s of 'in' and 'other' do not agree.",key);
			snprintf(key,3,"o%d",axis);
			if (!sf_histfloat(other,key,&oy2)) oy2=0.;
		}
		
		oy2 = -(oy2+(m2-1)*d2);
	} else {
		autocorr = true;
		one = false;
		m2 = n2;
		oy2 = -(o2+(n2-1)*d2);
	}


	
	nm2 = n2 + m2 - 1;
	snprintf(key,3,"n%d",axis);  sf_putint  (corr,key,nm2);
	snprintf(key,3,"o%d",axis);  sf_putfloat(corr,key,o2+oy2);
	
	// set up FFT 
	nfft  = 2*kiss_fft_next_fast_size((nm2+1)/2);
	nw    = nfft/2+1;
	scale = 1./nfft;
	cfg   = kiss_fftr_alloc(nfft,0,NULL,NULL);
	icfg  = kiss_fftr_alloc(nfft,1,NULL,NULL);
	fft1  = sf_complexalloc(nw);

	z  = sf_floatalloc(nfft*n1);
	z1 = sf_floatalloc(n1);
	
	n12 = n1*n2;
	x  = sf_floatalloc(n12);
	x1 = sf_floatalloc(nfft); 
	for (i2=n2; i2<nfft; i2++) 
		x1[i2] = 0.;

	if (!autocorr) {	// crosscorrelation: allocate y, y1, fft2
		fft2  = sf_complexalloc(nw);
		m12 = one? m2:n1*m2;
		y1 = sf_floatalloc(nfft); 
		if (!one) 
			y = sf_floatalloc(m12);
		for (i2=m2; i2<nfft; i2++) 
			y1[i2] = 0.;
	
		if (one) {
			sf_floatread(y1,m2,other);
			for (i2=0; i2<m2/2; i2++) {	// time-reverse y1
				tmp=y1[i2];  y1[i2]=y1[m2-1-i2]; y1[m2-i2-1]=tmp;
			}
			kiss_fftr(cfg, y1, (kiss_fft_cpx*) fft2);
		}
	} else {		// autocorrelation: allocate and make phase-shift array
		ce = sf_complexalloc(nw);
		for (iw=0; iw<nw; iw++) {
			shift = -2*SF_PI*iw*(n2-1)/nfft;
			ce[iw] = sf_cmplx(cosf(shift),sinf(shift));
		}
	}
	
	// main program
	if (autocorr) {
		for (i3=0; i3<n3; i3++) {
			sf_floatread(x, n12, in);
			
			for (i1=0; i1<n1; i1++) {
				for (i2=0; i2<n2; i2++)		// transpose
					x1[i2] = x[i1+i2*n1];
				kiss_fftr(cfg, x1, (kiss_fft_cpx*) fft1);
				
				
				for (iw=0; iw<nw; iw++) {
#ifdef SF_HAS_COMPLEX_H
					fft1[iw] = fft1[iw]*conjf(fft1[iw])*ce[iw]*scale;
					
#else
					fft1[iw] = sf_crmul(sf_cmul(sf_cmul(
									fft1[iw], sf_conjf(fft1[iw])), ce[iw]), scale);
#endif
				}
				kiss_fftri(icfg,(kiss_fft_cpx*) fft1, z+i1*nfft);
			}
			
			for (iw=0; iw<nm2; iw++) {	// transpose back
				for (i1=0; i1<n1; i1++)
					z1[i1] = z[iw+nfft*i1];
				sf_floatwrite(z1,n1,corr);			
			}
			
		}
	} else {	//cross-correlate
		for (i3=0; i3<n3; i3++) {
			sf_floatread(x, n12, in);
			
			if (!one) sf_floatread(y,m12,other);
			
			for (i1=0; i1<n1; i1++) {
				for (i2=0; i2<n2; i2++)		// transpose
					x1[i2] = x[i1+i2*n1];
				kiss_fftr(cfg, x1, (kiss_fft_cpx*) fft1);
				
				if (!one) {
					for (i2=0; i2<m2; i2++) // transpose and time-reverse
						y1[m2-i2-1] = y[i1+i2*n1];
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
			
			for (iw=0; iw<nm2; iw++) {	// transpose back
				for (i1=0; i1<n1; i1++)
					z1[i1] = z[iw+nfft*i1];
				sf_floatwrite(z1,n1,corr);			
			}
			
		}
	}
	
    exit(0);
}


