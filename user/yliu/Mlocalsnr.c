/* Local signal-to-noise ratio (SNR) estimation. */
/*
  Copyright (C) 2017 Jilin University
  
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

int main (int argc, char* argv[]) 
{
    bool verb,stack;
    int n1, n2, n3;
    int i,i1, i2, i3, dim2,m, k,l, xc,yc,nw,n12,niter,rect[2],n[2],fold;
    float *input, *output,*noise,*num,*den,*rat,temp,eps;
    sf_file in,en,out;
    char key[6];

    sf_init (argc, argv); 
    in = sf_input("in");
    en = sf_input("en");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(en,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(en,"n2",&n2)) sf_error("No n2= in input");

    n3 = sf_leftsize(en,2);

    if (!sf_getint("nw",&nw)) sf_error("nw needs integer input");
    /*window length*/
    if (nw%2 ==0) nw +=1;
    m = (nw-1)/2;
   
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of inversion iterations */
    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */
    if (!sf_getbool("stack",&stack)) stack = true;
    /* if y, window centre point, whereas window average*/

    input = sf_floatalloc(n1*n2);
    noise = sf_floatalloc(n1*n2);
    output = sf_floatalloc(n1*n2);
    num = sf_floatalloc(nw*nw);
    den = sf_floatalloc(nw*nw);
    rat = sf_floatalloc(nw*nw);

    dim2=2;
    n12 = 1;
    for (i=0; i < 2 ; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	n[i]=nw;
	n12 *= n[i];
    }
    
    sf_divn_init(dim2, n12, n, rect, niter, verb);
    
    for(i3=0; i3 < n3; i3++) {
	sf_floatread(input,n1*n2,in);
	sf_floatread(noise,n1*n2,en);
	for(i2=0; i2 < n2; i2++) {
	    for(i1=0; i1 < n1; i1++) {
		for(k=0;k<nw;k++){
		    for(l=0;l<nw;l++){
			yc=i2-m+k;
			xc=i1-m+l;
			if (xc < 0 || xc >= n1 || yc < 0 || yc >= n2) {
			    num[k*nw+l]=0;
			    den[k*nw+l]=0;
			} else {
			    num[k*nw+l]=input[yc*n1+xc]*input[yc*n1+xc];
			    den[k*nw+l]=noise[yc*n1+xc]*noise[yc*n1+xc];
			}
		    }
		}
		
		sf_divne (num,den,rat,eps);
		
		if (stack) {
		    output[i2*n1+i1] =10*log10(rat[nw*m+m]);
		} else {
		    temp=0.;
		    fold=0;
                    for(k=0; k < nw*nw; k++) {
			temp += rat[k];
			if (0!=rat[k]) fold++;
		    }
		    temp /= (fold+FLT_EPSILON);
		    output[i2*n1+i1] =10*log10(temp);
		}
	    }
	}
	sf_floatwrite(output,n1*n2,out);
    }
    
    exit(0);
}

/* 	$Id$	 */

