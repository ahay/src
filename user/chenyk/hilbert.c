/* Hilbert transform using different methods */
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
*/

#include<rsf.h>
#include"hilbert.h"
#include"fft.h"

void hilbert_t (float* trace /*input signal*/, 
		float* trace2 /*output hilbert data*/,
                int n1 /* trace length */)
/*< time domain tranform >*/
{
   int i,j;
   for(i=0;i<n1;i++)
	trace2[i]=0;
   for(i=0;i<n1;i++)
	for(j=0;j<n1-i;j++)
	     if(n1%2==1)
	    	 {trace2[i+j]+=trace[i]*(1.0/(tanf(SF_PI*j/n1)+SF_EPS)-cosf(SF_PI*j)/(sinf(SF_PI*j/n1)+SF_EPS))/n1; }
	     else
		 {trace2[i+j]+=trace[i]*sinf(SF_PI*j/2)*sinf(SF_PI*j/2)/(tanf(SF_PI*j/n1)+SF_EPS)*2/n1; }
}

void hilbert_f (float* trace /*input signal*/, 
		float* trace2 /*output hilbert data*/,
                int n1 /* trace length */)
/*< frequency domain transform >*/
{
    int j,nw;
    kiss_fft_scalar t;
    kiss_fft_cpx *ddf;
    nw=kiss_fft_next_fast_size((n1+1)/2)+1;
    ddf = (kiss_fft_cpx*) sf_complexalloc(nw);

    fft(trace,ddf,n1);
    for(j=0;j<nw;j++)
	{t=ddf[j].r; ddf[j].r=ddf[j].i; ddf[j].i=-t;}
    ifft(ddf,trace2,n1);
}

void hilbert_m (float* trace /*input signal*/, 
		float* trace2 /*output hilbert data*/,
                int nt1 /* trace length */, 
                int n1 /* transform length */, 
                int c1 /* filter parameter */)
/*< maximally flat hilbert transformer (MFHT) transform >*/
{
    int i,n,nt, it;
    float c,c2,*h;

    n = n1;
    nt = nt1;
    c = 1./(2*sqrtf(c1));
    c2 = c*c;
    h = sf_floatalloc(nt);

    for (it=0; it < nt; it++) {
	h[it] = trace[it];
    }

    for (i=n; i >= 1; i--) {
	trace2[0] = h[0] + (h[2]-2*h[1]+h[0])*c2;
	trace2[1] = h[1] + (h[3]-3*h[1]+2*h[0])*c2;
	for (it=2; it < nt-2; it++) {
	    trace2[it] = h[it]+(h[it+2]-2.*h[it]+h[it-2])*c2;
	}
	trace2[nt-2] = h[nt-2] + (h[nt-4]-3*h[nt-2]+2*h[nt-1])*c2;
	trace2[nt-1] = h[nt-1] + (h[nt-3]-2*h[nt-2]+h[nt-1])*c2;

	for (it=0; it < nt; it++) {
	    h[it] = trace[it] + trace2[it]*(2*i-1)/(2*i);
	}
    }

    trace2[0] = 2.*(h[0]-h[1])*c;
    for (it=1; it < nt-1; it++) {
	trace2[it] = (h[it-1]-h[it+1])*c;
    }

    trace2[nt-1] = 2.*(h[nt-2]-h[nt-1])*c;
}
