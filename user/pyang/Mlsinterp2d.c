/* Least-squares interpolation for 2D validition
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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

  Reference: On analysis-based two-step interpolation methods for randomly 
	sampled seismic data, P Yang, J Gao, W Chen, Computers & Geosciences
*/
#include <rsf.h>
#include <math.h>
#include <complex.h>

#include "wfftn.h"
#include "masklop.h"

int main(int argc, char* argv[])
{
    bool verb;
    int i, iter, niter, nouter, n[2]; 
    float eps, err;
    sf_file Fin, Fout, Fmask;/* mask and I/O files*/ 

    /* define temporary variables */
    int n1,n2;
    float *dat, *mask, *w;
    sf_complex *dobs, *dd, *p;

    sf_init(argc,argv);	/* Madagascar initialization */

    /* setup I/O files */
    Fin=sf_input("in");	/* read the data to be interpolated */
    Fmask=sf_input("mask");  	/* read the 2-D mask for 3-D data */
    Fout=sf_output("out"); 	/* output the reconstructed data */
 
    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity */
    if (!sf_getint("niter",&niter)) 	niter=100;
    if (!sf_getint("nouter",&nouter)) 	nouter=5;
    /* total number iterations */
    if (!sf_getfloat("eps",&eps)) 	eps=1.0e-2;
    /* regularization parameter */

    /* Read the data size */
    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");

    /* allocate data and mask arrays */
    dat=sf_floatalloc(n1*n2);
    dobs=sf_complexalloc(n1*n2);
    dd=sf_complexalloc(n1*n2);
    sf_floatread(dat,n1*n2,Fin);
    for(i=0; i<n1*n2; i++) dobs[i]=dat[i];
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2);
	sf_floatread(mask,n2,Fmask);
    } else {
	mask=NULL;
    }

    n[0]=n1; n[1]=n2;
    mask_init(n1, n2, mask);
    w=sf_floatalloc(n1*n2);

    for(i=0; i<n1*n2; i++) w[i]=10.0;
    fftn_init(2, n, w);
    if(1){
      for(iter=0; iter<nouter; iter++){
	/* each iteration reset w */
	if(iter==0) for(i=0; i<n1*n2; i++) w[i]=10.0;
	else {
	  fftn_fft(dd);
	  for(i=0; i<n1*n2; i++) w[i]=cabsf(dd[i]); 
	}
	/* each iteration reset dobs=dat and dd=0 */
	for(i=0; i<n1*n2; i++) dobs[i]=dat[i];
	for(i=0; i<n1*n2; i++) dd[i]=0.0;
	sf_csolver_prec(mask_lop, sf_ccgstep, fftn_lop, n1*n2, n1*n2, n1*n2, dd, dobs, niter, eps,"verb",verb,"end");
	/* for(i=0; i<n1*n2; i++) w[i]*=cabsf(p[i]); */
      }
    }

    fftn_close();
    sf_ccgstep_close();

    for(i=0; i<n1*n2; i++) dat[i]=crealf(dd[i]);
    sf_floatwrite(dat,n1*n2,Fout); /* output reconstructed seismograms */

    exit(0);
}
