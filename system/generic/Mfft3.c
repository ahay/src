/* FFT transform on extra axis.

Input and output are complex data. The input is padded by factor pad.
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "fftlabel.h"

int main (int argc, char **argv)
{
    int n1, nx, n3, dim, n[SF_MAX_DIM];  /* dimensions */
    int i1, ix, i3, j;       /* loop counters */
    int nk;                  /* number of wavenumbers */   
    int npad;                /* padding */

    float dx;                /* space sampling interval */
    float dk;                /* wavenumber sampling interval */
    float x0;                /* staring space */
    float k0;                /* starting wavenumber */
    float wt;                /* Fourier scaling */

    kiss_fft_cpx **cp, *ctrace; /* frequency-wavenumber */

    bool inv;                /* forward or inverse */
    bool sym;                /* symmetric scaling */
    bool opt;                /* optimal padding */
    int sign;                /* transform sign */
    int axis;                /* transform axis */

    char varname[12];        /* variable name */
    char *label;             /* transformed axis label */
    kiss_fft_cfg cfg;
    
    sf_file in, out;
    
    sf_init(argc,argv);
    in  = sf_input ( "in");
    out = sf_output("out");
    
    if (SF_COMPLEX != sf_gettype(in)) sf_error ("Need complex input");
    
    if (!sf_getbool("inv",&inv)) inv = false;
    /* if y, perform inverse transform */
    
    if (!sf_getbool("sym",&sym)) sym=false;
    /* if y, apply symmetric scaling to make the FFT operator Hermitian */
    
    if (!sf_getint("sign",&sign)) sign = inv? 1: 0;
    /* transform sign (0 or 1) */
    
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    
    if (!sf_getint("axis",&axis)) axis=2;
    /* Axis to transform */
    
    dim = sf_filedims(in,n);
    
    n1=n3=1;
    for (j=0; j < dim; j++) {
	if      (j < axis-1) n1 *= n[j];
	else if (j > axis-1) n3 *= n[j]; 
    }
    
    if (inv) { 
	sprintf(varname,"n%d",axis);
	if (!sf_histint  (in,varname,&nk)) sf_error("No %s= in input",varname);
	sprintf(varname,"d%d",axis);
	if (!sf_histfloat(in,varname,&dk)) sf_error("No %s= in input",varname);
	
	sprintf(varname,"fft3_n%d",axis);
	if (!sf_histint  (in,varname,&nx)) nx=nk;
	sprintf(varname,"fft3_o%d",axis);
	if (!sf_histfloat(in,varname,&x0)) x0 = 0.; 
	sprintf(varname,"fft3_label%d",axis);
	label = sf_histstring(in,varname);
	
	dx = 1./(nk*dk);
	
	sprintf(varname,"n%d",axis);
	sf_putint (out,varname,nx);
	sprintf(varname,"d%d",axis);
	sf_putfloat (out,varname,dx);
	sprintf(varname,"o%d",axis);
	sf_putfloat (out,varname,x0);
	sprintf(varname,"label%d",axis);
	if (NULL != label) {
	    sf_putstring(out,varname,label);
	} else if (NULL != (label = sf_histstring(in,varname))) {
	    (void) fix_label(axis,label,out);
	}
    } else { 
	sprintf(varname,"n%d",axis);
	if (!sf_histint  (in,varname,&nx)) sf_error("No %s= in input",varname);
	sprintf(varname,"d%d",axis);
	if (!sf_histfloat(in,varname,&dx)) sf_error("No %s= in input",varname);
	sprintf(varname,"o%d",axis);
	if (!sf_histfloat(in,varname,&x0)) x0 = 0.;
	sprintf(varname,"label%d",axis);
	label = sf_histstring(in,varname);
	
	sprintf(varname,"fft3_n%d",axis);
	sf_putint(out,varname,nx);
	sprintf(varname,"fft3_o%d",axis);
	sf_putfloat(out,varname,x0);
	if (NULL != label) {
	    sprintf(varname,"fft3_label%d",axis);
	    sf_putstring(out,varname,label);
	}
	
	if (!sf_getint("pad",&npad)) npad=2;
	/* padding factor */
	
	/* determine wavenumber sampling */
	nk = opt? kiss_fft_next_fast_size(nx*npad): nx*npad;
	if (nk != nx) sf_warning("padded to %d",nk);
	
	dk = 1./(nk*dx);
	k0 = -0.5/dx;
	
	sprintf(varname,"n%d",axis);
	sf_putint (out,varname,nk);
	sprintf(varname,"d%d",axis);
	sf_putfloat (out,varname,dk);
	sprintf(varname,"o%d",axis);
	sf_putfloat (out,varname,k0);
	if (NULL != label && !fix_label(axis,label,out)) {
	    sprintf(varname,"label%d",axis);
	    sf_putstring(out,varname,"Wavenumber");
	}
    }
    fix_unit(axis,in,out);
    
    cfg = kiss_fft_alloc(nk,sign,NULL,NULL);
    
    cp     = (kiss_fft_cpx**) sf_complexalloc2(n1,nk);
    ctrace = (kiss_fft_cpx*)  sf_complexalloc(nk);
    
    /* FFT scaling */
    wt = sym? 1./sqrtf((float) nk): 1./nk;
    
    for (i3=0; i3<n3; i3++) {
	if (inv) {
	    sf_floatread((float*) cp[0],n1*nk*2,in);
	    
	    for (i1=0; i1 < n1; i1++) {
		/* Fourier transform k to x */
		kiss_fft_stride(cfg,cp[0]+i1,ctrace,n1);
		
		for (ix=0; ix<nx; ix++) {
		    cp[ix][i1] = sf_crmul(ctrace[ix],ix%2? -wt: wt);
		}
	    }
	    
	    sf_floatwrite((float*) cp[0],n1*nx*2,out);
	} else {
	    sf_floatread((float*) cp[0],n1*nx*2,in);
	    
	    /* FFT centering */
	    for (ix=1; ix<nx; ix+=2) {
		for (i1=0; i1<n1; i1++) {
		    cp[ix][i1] = sf_cneg(cp[ix][i1]);
		}
	    }
	    
	    if (sym) {
		for (ix=0; ix<nx; ix++) {
		    for (i1=0; i1 < n1; i1++) {
			cp[ix][i1] = sf_crmul(cp[ix][i1],wt);
		    }
		}
	    }
	    
	    /* pad with zeros */
	    for (ix=nx; ix<nk; ix++) {
		for (i1=0; i1<n1; i1++) {
		    cp[ix][i1].r = 0.;
		    cp[ix][i1].i = 0.;
		}
	    }
	    
	    for (i1=0; i1 < n1; i1++) {
		/* Fourier transform x to k */
		kiss_fft_stride(cfg,cp[0]+i1,ctrace,n1);
		
		/* Transpose */
		for (ix=0; ix<nk; ix++) {
		    cp[ix][i1] = ctrace[ix];
		}
	    }
	    
	    sf_floatwrite((float*) cp[0],n1*nk*2,out);
	}
    }
    
    sf_fileclose(in);

    exit (0);
}

/* 	$Id$	 */
