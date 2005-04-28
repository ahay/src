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

int main (int argc, char **argv)
{
    int n1, nx, n3, axis, dim, n[SF_MAX_DIM];	/* dimensions */
    int i1, ix, i3, j;       /* loop counters 	*/
    int nk;		     /* number of wavenumbers */	
    int npad;                /* padding */

    float dx;		     /* space sampling interval */
    float dk;	             /* wavenumber sampling interval */
    float x0;                /* staring space */
    float k0;                /* starting wavenumber */
    float wt;                /* Fourier scaling */

    float complex **cp, *ctrace; /* frequency-wavenumber */

    bool inv;                /* forward or inverse */
    int sign;                /* transform sign */

    char varname[6];         /* variable name */
    kiss_fft_cfg cfg;
    
    sf_file in, out;

    sf_init(argc,argv);
    in  = sf_input ( "in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error ("Need complex input");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* if y, perform inverse transform */

    if (!sf_getint("sign",&sign)) sign = inv? 1: 0;
    /* transform sign (0 or 1) */

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

	sprintf(varname,"m%d",axis);
	if (!sf_histint  (in,varname,&nx)) nx=nk;
	sprintf(varname,"c%d",axis);
	if (!sf_histfloat(in,varname,&x0)) x0 = 0.; 

	dx = 1./(nk*dk);

	sprintf(varname,"n%d",axis);
	sf_putint (out,varname,nx);
	sprintf(varname,"d%d",axis);
	sf_putfloat (out,varname,dx);
	sprintf(varname,"o%d",axis);
	sf_putfloat (out,varname,x0);
    } else { 
	sprintf(varname,"n%d",axis);
	if (!sf_histint  (in,varname,&nx)) sf_error("No %s= in input",varname);
	sprintf(varname,"d%d",axis);
	if (!sf_histfloat(in,varname,&dx)) sf_error("No %s= in input",varname);
	sprintf(varname,"o%d",axis);
	if (!sf_histfloat(in,varname,&x0)) x0 = 0.;

	sprintf(varname,"m%d",axis);
	sf_putint(out,varname,nx);
	sprintf(varname,"c%d",axis);
	sf_putfloat(out,varname,x0);

	if (!sf_getint("pad",&npad)) npad=2;
	/* padding factor */

	/* determine wavenumber sampling */
	nk = nx*npad;
	dk = 1./(nk*dx);
	k0 = -0.5/dx;

	sprintf(varname,"n%d",axis);
	sf_putint (out,varname,nk);
	sprintf(varname,"d%d",axis);
	sf_putfloat (out,varname,dk);
	sprintf(varname,"o%d",axis);
	sf_putfloat (out,varname,k0);
    }
    
    cfg = kiss_fft_alloc(nk,sign,NULL,NULL);

    cp = sf_complexalloc2(n1,nk);
    ctrace = sf_complexalloc(nk);

    for (i3=0; i3<n3; i3++) {
	if (inv) {
	    sf_complexread(cp[0],n1*nk,in);

	    /* FFT scaling */
	    wt = 1./nk;

	    for (i1=0; i1 < n1; i1++) {
		/* Fourier transform k to x */
		kiss_fft_stride(cfg,(kiss_fft_cpx *) (cp[0]+i1),
				(kiss_fft_cpx *) ctrace,n1);
		
		for (ix=0; ix<nx; ix++) {
		    cp[ix][i1] = ix%2? -ctrace[ix]*wt: ctrace[ix]*wt;
		}
	    }
      
	    sf_complexwrite(cp[0],n1*nx,out);
	} else {
	    sf_complexread(cp[0],n1*nx,in);
	    
	    /* FFT centering */
	    for (ix=1; ix<nx; ix+=2) {
		for (i1=0; i1<n1; i1++) {
		    cp[ix][i1] = -cp[ix][i1];
		}
	    }

	    /* pad with zeros */
	    for (ix=nx; ix<nk; ix++) {
		for (i1=0; i1<n1; i1++) {
		    cp[ix][i1] = 0.;
		}
	    }
    
	    for (i1=0; i1 < n1; i1++) {
		/* Fourier transform x to k */
		kiss_fft_stride(cfg,(kiss_fft_cpx *) (cp[0]+i1),
				(kiss_fft_cpx *) ctrace,n1);
		/* Transpose */
		for (ix=0; ix<nk; ix++) {
		    cp[ix][i1] = ctrace[ix];
		}
	    }
	    
	    sf_complexwrite(cp[0],n1*nk,out);
	}
    }
    
    exit (0);
}

/* 	$Id$	 */
