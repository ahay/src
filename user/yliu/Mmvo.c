/* Calculate MVO and PVO curve of CSEM data. */
/*
  Copyright (C) 2013 Jilin University
  
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
    bool opt, mvo, log;
    int n1, nw1, nt, nw, i, j, N, k, n, nnw, m1;
    float d1, f, dw;
    float *um, *TD, *outp;
    kiss_fft_cpx *FD;
    kiss_fftr_cfg cfg;
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
    nt = opt? 2*kiss_fft_next_fast_size((nw+1)/2): nw;

    if (nt%2) nt++;
    nw1 = nt/2+1;
    dw = 1./(nt*d1);
    k=(int)(f/dw);
    if(fabs(f-k*dw)>fabs(f-(k+1)*dw))
        k=k+1;

    FD = (kiss_fft_cpx*)sf_complexalloc(nw1);
    TD = sf_floatalloc(nt);
    cfg = kiss_fftr_alloc(nt,0,NULL,NULL);

    for(m1=n1; ;m1++) {
	if(m1%nw==0) break;
    }
    N=m1/nw;

    sf_putint(out, "n1", N);
    sf_putfloat(out, "o1", 0);
    sf_putfloat(out, "d1", nw*d1);

    outp =sf_floatalloc(N);

    for(j=0;j<N;j++) {
	/*  sf_warning("slice %d of %d;",j+1,N); */
	
	for(i=0;i<nt;i++){
	    if((i<nw)&&((j*nw+i)<n1)) {
		TD[i]=um[j*nw+i];
	    }
	    else{TD[i]=0.;} 
	}
	
	kiss_fftr(cfg,TD,FD);  
	
	if(mvo) {
	    if(log) {
		outp[j]=log10(sf_cabsf(FD[k])*2/nw);
	    } else {
		outp[j]=sf_cabsf(FD[k])*2/nw;
	    }
	}
	if(!mvo) outp[j]= 180.0*sf_cargf(FD[k])/SF_PI;
    }    

    sf_floatwrite(outp,N,out);
    
    free(um);
    free(TD);
    free(outp);
    
    exit(0);
}

/* 	$Id$	 */
