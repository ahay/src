/* Calculate MVO and PVO curve of CSEM data (another version). */
/*
  Copyright (C) 2014 Jilin University
  
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
#include <stdio.h>
#include <math.h>

void FFT(float x[], float y[], int n, int sign)
{
    int i,j,k,l,m,n1,n2;
    float c,c1,e,s,s1,t,tr,ti;
    
    /* Calculate i = log2N */
    for(j = 1,i = 1; i<1000; i++) {
	m = i;
	j = 2*j;
	if(j == n)
	    break;
    }
    
    n1 = n - 1;
    for(j=0,i=0; i<n1; i++) {
	if(i<j) {
	    tr = x[j];
	    ti = y[j];
	    x[j] = x[i];
	    y[j] = y[i]; 
	    x[i] = tr;
	    y[i] = ti;                 
	}
	k = n/2;
	while(k<(j+1)) {
	    j = j - k;
	    k = k/2;              
	}
	j = j + k;
    }
    
    n1 = 1;
    for(l=1; l<=m; l++) {
	n1 = 2*n1;
	n2 = n1/2;
	e = SF_PI/n2;
	c = 1.0;
	s = 0.0;
	c1 = cos(e);
	s1 = -sign*sin(e);
	for(j=0; j<n2; j++) {
	    for(i=j; i<n; i+=n1) {
		k = i + n2;
		tr = c*x[k] - s*y[k];
		ti = c*y[k] + s*x[k];
		x[k] = x[i] - tr;
		y[k] = y[i] - ti;
		x[i] = x[i] + tr;
		y[i] = y[i] + ti;        
	    }
	    t = c;
	    c = c*c1 - s*s1;
	    s = t*s1 + s*c1;
	}
    }
    if(sign == -1) {
	for(i=0; i<n; i++) {
	    x[i] /= n;
	    y[i] /= n;
	}
    }
}

int main (int argc, char* argv[])
{
    bool opt, mvo, log;
    int n1, nt, nw, i, j, N, k, n, nnw, m1;
    float d1, f, dw;
    float *um, *TD, *outp;
    int fftk;
    int fftn;
    float *FD;
    sf_file in, out;

    sf_init (argc, argv);

    in= sf_input("in");
    out = sf_output("out");
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */

    if (!sf_getbool("log",&log)) log=true;
    /* if y, calculate logarithm of MVO */

    if (!sf_getbool("mvo",&mvo)) mvo=true;
    /* if y, MVO curve; otherwise, PVO curve */    

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_getfloat("f",&f)) f=0.08;
    /* calculate frequency */

    if (!sf_getint("nnw",&nnw)) nnw=1600;
    /* sample window */

    if (!sf_getint("n",&n)) n=1;
    /* number of window period */

    um = sf_floatalloc(n1);	
    nw=nnw*n;
    sf_floatread(um,n1,in);     

    
    fftk=(int)(log10((float)(nw))/log10((float)(2)));
    if(nw==pow(2,fftk))	{
	fftn=pow(2,fftk);
    } else {
	fftk++;
	fftn=pow(2,fftk);
    }
    
    nt = fftn; 

    if (nt%2) nt++;
    dw = 1./(nt*d1);
    k=(int)(f/dw);

    TD = sf_floatalloc(nt);
    FD = sf_floatalloc(nt);

  
    for(m1=n1; ;m1++) {
	if(m1%nw==0) break;
    }
    N=m1/nw;

    sf_putint(out, "n1", N);
    sf_putfloat(out, "o1", 0);
    sf_putfloat(out, "d1", 1);

    outp =sf_floatalloc(N);

    for(j=0;j<N;j++) {
	/*  sf_warning("slice %d of %d;",j+1,N); */
	
	for(i=0;i<nt;i++){
	    if((i<nw)&&((j*nw+i)<n1)) {
		TD[i]=um[j*nw+i];
		FD[i]=0.;
	    }
	    else{
		TD[i]=0.;
		FD[i]=0.;
	    } 
	}
	
	FFT(TD, FD, nt, 1);
	
	if(mvo) {
	    if(log) {
		outp[j]=log10(sqrt(TD[k]*TD[k]+FD[k]*FD[k]));
	    } else {
		outp[j]=sqrt(TD[k]*TD[k]+FD[k]*FD[k]);
	    }
	}
	if(!mvo) outp[j]=180.0*atan2(FD[k],TD[k])/SF_PI;
    }
    
    sf_floatwrite(outp,N,out);
    
    free(um);
    free(TD);
    free(outp);
    
    exit(0);
}

/* 	$Id$	 */
