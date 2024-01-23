/* Estimate local centroid frequency. */
/*
  Copyright (C) 2022 Jilin University
  
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

int main(int argc, char* argv[])
{
    bool  verb;

    int   n1, n2;                    
    float d1, o1;
    int   i, dim, nd, i1, i2, iw, nw, n1w, niter;
    int   n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float w, w0, dw, eps;
    float *ltftm, *inf, *infnum, *infden, *infvar2, *infvar;
    float *varn, *vard, *bw1, *bw2, *tran, *bran;
    sf_complex *ltft,*tfl,*tfh,*tfm;
    sf_file in,out,avef,var2,tflo,tfhi,tfme,trange,brange;
    char key[6];

    sf_axis ag1, ag2;
   
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    /* verbosity flag */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* smoothing radius */
    if (!sf_getint("niter",&niter)) niter = 100;
    /* regularization */
    if (!sf_getfloat("eps",&eps)) eps=0.0f;

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1. ;
    if (!sf_histfloat(in,"o1",&o1)) o1 = 0. ;

    n2=sf_leftsize(in,2);
    if(!sf_histint(in,"n2",&nw)) sf_error("No n2=in input");
    if(!sf_histfloat(in,"d2",&dw)) sf_error("No d2=in input");
    if(!sf_histfloat(in,"o2",&w0)) sf_error("No o2=in input");
    /* sf_unshiftdim(in, out, 1); */
    
    sf_warning("n2=%d",n2);
    
    if (NULL != sf_getstring("trange")) { /* data weight */
	trange = sf_input("trange");
	tran = sf_floatalloc(n1);	
    } else {
	trange = NULL;
	tran = NULL;
    }

    if (NULL != sf_getstring("brange")) { /* data weight */
	brange = sf_input("brange");
	bran = sf_floatalloc(n1);	
    } else {
	brange = NULL;
	bran = NULL;
    }
    
    if(NULL!=sf_getstring("avef")) {
    	avef = sf_output("avef");
    	sf_settype(avef,SF_FLOAT);
	
    } else {
    	avef = NULL;}

    if(NULL!=sf_getstring("var2")) {
    	var2 = sf_output("var2");
    	sf_settype(var2,SF_FLOAT);
	
    } else {
    	var2 = NULL;}
    
    if(NULL!=sf_getstring("tflo")) {
    	tflo = sf_output("tflo");
    	sf_settype(tflo,SF_COMPLEX);
    } else {
    	tflo = NULL;}

   if(NULL!=sf_getstring("tfhi")) {
    	tfhi = sf_output("tfhi");
    	sf_settype(tfhi,SF_COMPLEX);
    } else {
    	tfhi = NULL;}

   if(NULL!=sf_getstring("tfme")) {
    	tfme = sf_output("tfme");
    	sf_settype(tfme,SF_COMPLEX);
    } else {
    	tfme = NULL;}
    
    /* set axis */
    ag1 = sf_maxa(n1,o1,d1);
    ag2 = sf_maxa(2, 0., 1);
    
    
    sf_oaxa(out,ag1,1);
    sf_oaxa(out,ag2,2);
    

    sf_settype(out,SF_FLOAT);

    if(NULL!=avef) {
    	sf_putint(avef,"n1",n1);
    	sf_putfloat(avef,"d1",d1);
    	sf_putfloat(avef,"o1",o1);
	sf_unshiftdim(in, avef, 2);
    }

    if(NULL!=var2) {
    	sf_putint(var2,"n1",n1);
    	sf_putfloat(var2,"d1",d1);
    	sf_putfloat(var2,"o1",o1);
	sf_unshiftdim(in, var2, 2);
    }
    
    n1w=n1*nw;    

    inf=sf_floatalloc(n1);
    infnum=sf_floatalloc(n1);
    infden=sf_floatalloc(n1);
    infvar2=sf_floatalloc(n1);
    infvar=sf_floatalloc(n1);
    varn=sf_floatalloc(n1);
    vard=sf_floatalloc(n1);
    bw1=sf_floatalloc(n1);
    bw2=sf_floatalloc(n1);
 
    ltft=sf_complexalloc(n1w);
    ltftm=sf_floatalloc(n1w);
    tfl = sf_complexalloc(n1w);
    tfm = sf_complexalloc(n1w);
    tfh = sf_complexalloc(n1w);
    
    /* bandwidth and local frequency */
    
    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);

	if(NULL != trange) {
	    sf_floatread(tran,n1,trange);
	}
	if(NULL != brange) {
	    sf_floatread(bran,n1,brange);
	}
	
	for (i1=0; i1 < n1; i1++) {
	    for (iw=0; iw < nw; iw++) {
		ltftm[i1*nw+iw]=0.;
	    }
	}
	sf_complexread(ltft,n1w,in);
	
	for (i1=0; i1 < n1; i1++) {
	    for (iw=0; iw < nw; iw++) {
		ltftm[i1*nw+iw]=sqrtf(crealf(ltft[iw*n1+i1])*
				      crealf(ltft[iw*n1+i1])
				      +cimagf(ltft[iw*n1+i1])*
				      cimagf(ltft[iw*n1+i1]));
	    }
	}

	dim = 1;
	nd =1;
	for (i=0; i < dim ; i++) {
	    snprintf(key,6,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=10;
            /*( rect#=(10,10,...) smoothing radius on #-th axis )*/	    
	    n[i] = n1;
	    nd *=n[i];	    
	}
	sf_divn_init(dim, nd, n, rect, niter, verb);
	sf_divn_init(dim, nd, n, rect, niter, verb);
	
	/* centroid frequency */
	for (i1=0; i1 < n1; i1++) {
	    inf[i1]=0.;
	    infnum[i1]=0.;
	    infden[i1]=0.;
	    for (iw=0; iw < nw; iw++) {
		w = w0 + iw*dw;
		if (NULL != trange && NULL == brange) {
		    if (w <= tran[i1]) {
			infnum[i1]+=w*ltftm[i1*nw+iw];
			infden[i1]+=ltftm[i1*nw+iw];			
		    }
		} else if (NULL != trange && NULL != brange) {
		    if (w <= tran[i1] && w >= bran[i1]) {
			infnum[i1]+=w*ltftm[i1*nw+iw];
			infden[i1]+=ltftm[i1*nw+iw];			
		    }
		} else {
		    infnum[i1]+=w*ltftm[i1*nw+iw];
		    infden[i1]+=ltftm[i1*nw+iw];
		}
	    }
	}	
	
	sf_divne(infnum, infden, inf, eps);

	if(NULL!=avef) sf_floatwrite(inf,n1,avef);

	/* spectral bandwidth */
	
	for (i1=0; i1 < n1; i1++) {
	    varn[i1]=0.;
	    vard[i1]=0.;
	    infvar2[i1]=0.;
	    infvar[i1]=0.;
	    for (iw=0; iw < nw; iw++) {
		w = w0 + iw*dw;
		if (NULL != trange && NULL == brange) {
		    if (w <= tran[i1]) {			
			varn[i1]+=(w-inf[i1])*(w-inf[i1])*ltftm[i1*nw+iw];
			vard[i1]+=ltftm[i1*nw+iw];
		    }
		} else if (NULL != trange && NULL != brange) {
		    if (w <= tran[i1] && w >= bran[i1]) {
			varn[i1]+=(w-inf[i1])*(w-inf[i1])*ltftm[i1*nw+iw];
			vard[i1]+=ltftm[i1*nw+iw];			
		    }
		} else {
		    varn[i1]+=(w-inf[i1])*(w-inf[i1])*ltftm[i1*nw+iw];
		    vard[i1]+=ltftm[i1*nw+iw];
		}		
	    }
	}	
		
	sf_divne(varn, vard, infvar2, eps);
	
	if(NULL!=var2) sf_floatwrite(infvar2,n1,var2);
	
	for (i1=0; i1 < n1; i1++) {
	    infvar[i1]=sqrt(infvar2[i1]);
	}
	
	for (i1=0; i1 < n1; i1++) {
	    bw1[i1]=inf[i1]-infvar[i1];
	    bw2[i1]=inf[i1]+infvar[i1];
	    if(bw1[i1]<0) {
		bw1[i1]=0.0;}
	    if(bw2[i1]<0) {
		bw2[i1]=0.0;}
	}
	sf_floatwrite(bw1,n1,out);
	sf_floatwrite(bw2,n1,out);
	
        /* spectral decomposition */
	for (iw=0; iw < nw; iw++) {
	    for (i1=0; i1<n1; i1++) {
		w = w0 + iw*dw;
		if(NULL!=tflo) {
		    if(w<bw1[i1]) {
			tfl[iw*n1+i1]=ltft[iw*n1+i1];
		    } else {
			tfl[iw*n1+i1]=0;}
		} 

		if(NULL!=tfhi) {
		    if(w>bw2[i1]) {
			tfh[iw*n1+i1]=ltft[iw*n1+i1];
		    } else {
			tfh[iw*n1+i1]=0;}
		} 

		if(NULL!=tfme) {
		    if(w>=bw1[i1] && w<=bw2[i1]) {
			tfm[iw*n1+i1]=ltft[iw*n1+i1];
		    } else {
			tfm[iw*n1+i1]=0;}
		}
	    }
	}
	if(NULL!=tflo) sf_complexwrite(tfl,n1w,tflo);
	if(NULL!=tfhi) sf_complexwrite(tfh,n1w,tfhi);
	if(NULL!=tfme) sf_complexwrite(tfm,n1w,tfme);
    }
    sf_warning(".");

    exit(0);
}
