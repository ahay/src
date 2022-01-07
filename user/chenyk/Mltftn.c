/* Non-stationary local time-frequency transform (NLTFT). 
*/
/*
  Copyright (C) 2021 University of Texas at Austin
  
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
#include "divnn_sc.h"

int main(int argc, char* argv[])
{
    bool inv, verb;
    int i1, n1, iw, nt, nw, i2, n2, rect0, rectf, niter, n12, n1w;
    int m[SF_MAX_DIM], *rect;
    float t, d1, w, w0, dw, mean=0.0f, alpha;
    float *trace, *kbsc, *mkbsc, *sscc=NULL, *mm, *ww;
    sf_complex *outp, *cbsc;
    sf_file in, out, basis;
   
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (NULL != sf_getstring("basis")) {
	basis = sf_output("basis");
	sf_settype(basis,SF_COMPLEX);
    } else {
	basis = NULL;
    }

    if (!inv) {
	if (!sf_getint("nw",&nw)) { /* number of frequencies */
	    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
	    nw = nt/2+1;
	    dw = 1./(nt*d1);
	    w0 = 0.;
	} else {
	    if (!sf_getfloat("dw",&dw)) {
		/* frequency step */
		nt = 2*kiss_fft_next_fast_size((n1+1)/2);
		dw = 1./(nt*d1);
	    }
	    if (!sf_getfloat("w0",&w0)) w0=0.;
	    /* first frequency */
	}
	n2 = sf_leftsize(in,1);


	if (!sf_getint("rect",&rect0)) rect0=10;
	/* smoothing radius (in time, samples) */
	if (!sf_getint("rectf",&rectf)) rectf=5;
	/* smoothing radius (in frequency, samples) */
	
	if (!sf_getint("niter",&niter)) niter=100;
	/* number of inversion iterations */
	if (!sf_getfloat("alpha",&alpha)) alpha=0.;
	/* frequency adaptivity */

	sf_warning("Check 000");
	
	for(i2=0; i2 < SF_MAX_DIM; i2 ++) {
	    m[i2] = 1;
	}
	m[0] = n1;
	m[1] = nw;
	m[2] = n2;
    } else {
	n2 = sf_leftsize(in,2);
	if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&w0)) sf_error("No o2= in input");
	sf_unshiftdim(in, out, 2);
	sf_settype(out,SF_FLOAT);
    }

    if (NULL != basis) {
	sf_shiftdim(in, basis, 2);
	sf_putint(basis,"n2",nw);
	sf_putfloat(basis,"d2",dw);
	sf_putfloat(basis,"o2",w0);
	sf_putstring(basis,"label2","Frequency");
	sf_putstring(basis,"unit2","Hz");
    }

    n1w = n1*nw;
    n12 = 2*n1w;
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;

    rect = sf_intalloc(2*nw);
    for (iw=0; iw < nw; iw++) {
	rect[iw+nw] = rect[iw] = SF_MAX(1, (int) rect0/(1.0+alpha*iw/nw));
    }

	int i,ii, dim, n[SF_MAX_DIM], dim1, nd, b;
    int box[SF_MAX_DIM];
    char key[8];
    
	/*box means maximum (edge) padding for triangle smoothing, box[i]=max(rect[i])*/
	sf_file rects[SF_MAX_DIM], shift[SF_MAX_DIM]; 
    int   *sft[SF_MAX_DIM];	/* storing non-stationary shifting size */
    float *rct[SF_MAX_DIM]; /* storing non-stationary smoothing radii */


    sf_warning("Check 0");


    dim = sf_filedims (in,n);
    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
    }
    
	/*Calculate dim1*/
    dim1 = -1;
    for (i=0; i < dim+1; i++) {
	snprintf(key,6,"rect%d",i); /*note I change here*/
	if (NULL != sf_getstring(key)) {
	    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ ) rect0 means frequency rect */
	    rects[i] = sf_input(key);
	    if (SF_FLOAT != sf_gettype(rects[i])) sf_error("Need float %s",key);
	    dim1 = i;
	    snprintf(key,8,"shift%d",i); /*note I change here*/
	    if (NULL != sf_getstring(key)) {
		/*( shift# shifting of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
		shift[i] = sf_input(key);
		if (SF_INT != sf_gettype(shift[i])) sf_error("Need int %s",key);
	    } else {
		shift[i] = NULL;
	    }
	} else {
	    rects[i] = NULL;
	    shift[i] = NULL;
	}
    }
	int n21;
	/*Calculate n1*/
	/*dim1: smoothing to dimension dim1*/
    n21 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i < dim1) {
	    n21 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

	n1w=n21*nw;
	n12=2*n1w;
    trace = sf_floatalloc(n21);
    kbsc    = sf_floatalloc(n12);
    outp = sf_complexalloc(n1w);
    cbsc = sf_complexalloc(n1w);

    sf_warning("n21=%d,n2=%d,n12=%d,n1=%d,dim=%d,dim1=%d",n21,n2,n12,n1,dim,dim1);
    sf_warning("m[0]=%d,m[1]=%d,m[2]=%d",m[0],m[1],m[2]);

	/*reading the non-stationary smoothing radii*/
    for (i=0; i < dim+1; i++) {
	if (NULL != rects[i]) {
	  if(i==0) {ii=dim1;}
	  else
	  {	ii=i-1;
// 	   if(i==1){ ii=0;}
// 	  	else
// 	  	{ii=i;}
	  }
		box[ii] = 1;
	  	rct[ii] = sf_floatalloc (n21*nw*2);
	    sft[ii] = sf_intalloc (n21*nw*2);

	    sf_floatread(rct[ii],n21*nw*2,rects[i]);
	    sf_fileclose(rects[i]);

	    if (NULL != shift[i]) {
		sf_intread(sft[ii],n21*nw*2,shift[i]);
		sf_fileclose(shift[i]);
	    } else {
		for (i1=0; i1 < n21*nw*2; i1++) {
		    sft[ii][i1] = 0;
		}
	    }

	    for (i1=0; i1 < n21*nw*2; i1++) {
		b = ceilf(rct[ii][i1])+SF_ABS(sft[ii][i1]);
		if (b > box[ii]) box[ii] = b;
	    }	   
	} else {
	    rct[i] = NULL;
	    sft[i] = NULL;
	}
    }

    int ntmp=0;
    ntmp=ceilf(n21/n1);
    
    sf_warning("box[0]=%d,box[1]=%d,box[2]=%d,box[3]=%d",box[0],box[1],box[2],box[3]);
    sf_warning("rct[0][0]=%g",rct[0][0]);
    sf_warning("rct[1][0]=%g",rct[1][0]);
    
    if(dim1>=2)
    sf_warning("rct[2][0]=%g",rct[2][0]);

    if(dim1>=3)
    sf_warning("rct[3][0]=%g",rct[3][0]);
    
    if (!inv) {

	/*dimension: n1*nw*2*/
	if(dim1==1)
	{
	m[0]=n1;
	m[1]=nw;
	m[2]=2;
	sscc = sf_floatalloc(n1*nw*2);
	}else{
	if(dim1==2)
	{m[0]=n1;
	m[1]=ntmp;
	m[2]=nw;
	m[3]=2;
	sscc = sf_floatalloc(n1*ntmp*nw*2);
	}else{
	if(dim1==3)
	{
	m[0]=n1;
	m[1]=n[1];
	m[2]=n[2];
	m[3]=nw;
	m[4]=2;
	sscc = sf_floatalloc(n1*ntmp*nw*2);
	}
	}
	}
	box[dim1+1]=1;
	divnn_sc_init2(2*nw, dim1+2, n21, m, box, rct,sft, kbsc, 
			(bool) (verb && (n2 < 500))); 
    } else {
	sscc = NULL;
    }
    
    sf_warning("m[0]=%d,m[1]=%d,m[2]=%d,m[3]=%d",m[0],m[1],m[2],m[3]);
    
    
    if (NULL != sf_getstring("mask")) { /* data weight */
	//mask = sf_input("mask");
	mm = sf_floatalloc(n1);	
    } else {
	//mask = NULL;
	mm = NULL;
    }

    if (NULL != sf_getstring("weight")) { /* model weight */
	//weight = sf_input("weight");
	ww = sf_floatalloc(n1w);
    } else {
	//weight = NULL;
	ww = NULL;
    }


	sf_shiftdim(in, out, 2);
	if(dim1==1)
	{
	sf_putint(out,"n2",nw);
	sf_putfloat(out,"d2",dw/(2.*SF_PI));
	sf_putfloat(out,"o2",w0/(2.*SF_PI));
	sf_putstring(out,"label2","Frequency");
	sf_putstring(out,"unit2","Hz");
	sf_settype(out,SF_COMPLEX);
	}
	else
	{
	if(dim1==2)
		{
		sf_putint(out,"n3",nw);
		sf_putfloat(out,"d3",dw/(2.*SF_PI));
		sf_putfloat(out,"o3",w0/(2.*SF_PI));
		sf_putstring(out,"label3","Frequency");
		sf_putstring(out,"unit3","Hz");
		sf_putint(out,"n2",ntmp);
		sf_putfloat(out,"d2",1);
		sf_putfloat(out,"o2",1);
		sf_settype(out,SF_COMPLEX);
		}else{
	if(dim1==3)
	{
		sf_putint(out,"n4",nw);
		sf_putfloat(out,"d4",dw/(2.*SF_PI));
		sf_putfloat(out,"o4",w0/(2.*SF_PI));
		sf_putstring(out,"label4","Frequency");
		sf_putstring(out,"unit4","Hz");
		sf_putint(out,"n3",n[2]);
		sf_putint(out,"n2",n[1]);
		sf_putfloat(out,"d2",1);
		sf_putfloat(out,"o2",1);
		sf_settype(out,SF_COMPLEX);
	
	}
	}
	}


    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
    for (i2=0; i2 < ntmp; i2++){
	for (i1=0; i1 < n1; i1++) {
	    if (0.==w) { /* zero frequency */
		kbsc[iw*ntmp*n1+i2*n1+i1] = 0.;
	    } else {
		t = i1*d1;
		kbsc[iw*ntmp*n1+i2*n1+i1] = sinf(w*t);
	    }
	}
    }
    }
    

    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
    for (i2=0; i2 < ntmp; i2++){
	for (i1=0; i1 < n1; i1++) {
	    if (0.==w) { /* zero frequency */
		kbsc[iw*ntmp*n1+i2*n1+i1+nw*ntmp*n1] = 0.5;
	    } else {
		t = i1*d1;
		kbsc[iw*ntmp*n1+i2*n1+i1+nw*ntmp*n1] = cosf(w*t);
	    }

	    cbsc[iw*ntmp*n1+i2*n1+i1] = sf_cmplx(kbsc[iw*ntmp*n1+i2*n1+i1+nw*ntmp*n1],
				      kbsc[iw*ntmp*n1+i2*n1+i1]);
	}
    }
	}
	sf_warning("ntmp=%d,nw=%d",ntmp,nw);
	
    if (NULL != mm || NULL != ww) {
	mkbsc = sf_floatalloc(n12);
	for (i1=0; i1 < n12; i1++) {
	    mkbsc[i1] = kbsc[i1];
	}
    } else {
	mkbsc = NULL;

	/*n12=N1*N2*NW*2=65016*/
	mean = 0.;
	for (i1=0; i1 < n12; i1++) {
	    mean += kbsc[i1]*kbsc[i1];
	}
	mean = sqrtf (mean/(n12));
	for (i1=0; i1 < n12; i1++) {
	    kbsc[i1] /= mean;
	}
    }

    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);

	if (NULL != basis) sf_complexwrite(cbsc,n1w,basis);

// 	if (NULL != mm || NULL != ww) {
// 	    for (i1=0; i1 < n12; i1++) {
// 		kbsc[i1] = mkbsc[i1];
// 	    }
// 
// 	    if (NULL != mm) {
// 		sf_floatread(mm,n21,mask);
// 		for (iw=0; iw < 2*nw; iw++) {
// 		    for (i1=0; i1 < n21; i1++) {
// 			kbsc[iw*n1+i1] *= mm[i1];
// 		    }
// 		}
// 	    }
// 
// 	    if (NULL != ww) {
// 		sf_floatread(ww,n1w,weight);
// 		for (iw=0; iw < nw; iw++) {
// 		    for (i1=0; i1 < n21; i1++) {
// 			kbsc[iw*n1+i1]      *= ww[iw*n1+i1];
// 			kbsc[(iw+nw)*n1+i1] *= ww[iw*n1+i1];
// 		    }
// 		}
// 	    }
// 
// 	    mean = 0.;
// 	    for (i1=0; i1 < n12; i1++) {
// 		mean += kbsc[i1]*kbsc[i1];
// 	    }
// 	    mean = sqrtf (mean/(n12));
// 	    for (i1=0; i1 < n12; i1++) {
// 		kbsc[i1] /= mean;
// 	    }
// 	}

	if (!inv) {
	    sf_floatread(trace,n21,in);
// 	    if (NULL != mm) {
// 		for (i1=0; i1 < n21; i1++) {
// 		    trace[i1] *= mm[i1];
// 		}
// 	    }
	    
	    for(i1=0; i1 < n21; i1++) {
		trace[i1] /= mean;
	    }
	    
	    /*Main program*/
	    divnn_sc2 (trace,sscc,niter);    
	    
	    for (iw=0; iw < nw; iw++) {
		for (i1=0; i1 < n21; i1++) {
		    outp[iw*n21+i1] = sf_cmplx(sscc[iw*n21+i1+nw*n21],
					      sscc[iw*n21+i1]);
		}
	    }

// 	    if (NULL != ww) {
// 		for (i1=0; i1 < n1w; i1++) {
// #ifdef SF_HAS_COMPLEX_H
// 		    outp[i1] *= ww[i1];
// #else
// 		    outp[i1] = sf_crmul(outp[i1],ww[i1]);
// #endif
// 		}
// 	    } 
		sf_warning("inversion done");
	    sf_complexwrite(outp,n1w,out);
	} else {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = 0.;
	    }
	    sf_complexread(outp,n1w,in);
	    for (iw=0; iw < nw; iw++) {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] += crealf(outp[iw*n1+i1])*kbsc[(iw+nw)*n1+i1]
			*mean+cimagf(outp[iw*n1+i1])*kbsc[iw*n1+i1]*mean;
		    if (NULL != mm) trace[i1] *= mm[i1];
		}
	    }
	    sf_floatwrite(trace,n1,out);
	}
    }
    sf_warning(".");

    exit(0);
}


