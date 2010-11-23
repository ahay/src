/* Multi-dimensional missing data interpolation. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

static sf_filter aa;
static float *temp;

static void polydiv_lop( bool adj, bool add, 
			 int nx, int ny, float* xx, float*yy) 
/* helical inverse convolution (polynomial division) */
{
    int ia, iy, ix;
    
    sf_adjnull( adj, add, nx, ny, xx, yy);
    
    for (ix=0; ix < nx; ix++) temp[ix] = 0.;
    
    if (adj) {
	for (ix = nx-1; ix >= 0; ix--) {  
	    temp[ix] = yy[ix];
	    for (ia = 0; ia < aa->nh; ia++) {
		iy = ix + aa->lag[ia];
		if( iy >= ny) continue;
		temp[ix] -= aa->flt[ia] * temp[iy];
	    }
	}
	for (ix=0; ix < nx; ix++) xx[ix] += temp[ix];
    } else {
	for (iy = 0; iy < ny; iy++) { 
	    temp[iy] = xx[iy];
	    for (ia = 0; ia < aa->nh; ia++) {
		ix = iy - aa->lag[ia];
		if( ix < 0) continue;
		temp[iy] -= aa->flt[ia] * temp[ix];
	    }
	}
	for (iy=0; iy < ny; iy++) yy[iy] += temp[iy];
    }
}

int main(int argc, char* argv[])
{
    int na,ia, niter, n,i;
    float a0, *mm, **xmov;
    bool prec, *known;
    char *lagfile;
    sf_file in, out, filt, lag, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    /* filter for inverse model covariance */
    out = sf_output("out");

    n = sf_filesize(in);

    if (!sf_getbool("prec",&prec)) prec=true;
    /* If y, use preconditioning */
    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1=");
    aa = sf_allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
 
    /* Get filter lags */
    if (NULL == (lagfile = sf_histstring(filt,"lag"))) {
	for (ia=0; ia < na; ia++) {
	    aa->lag[ia]=ia+1;
	}
    } else {
	lag = sf_input(lagfile);
	sf_intread(aa->lag,na,lag);
	sf_fileclose(lag);
    }

    /* Get filter values */
    sf_floatread(aa->flt,na,filt);
    sf_fileclose(filt);

    /* Normalize */
    for (ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    mm = sf_floatalloc(n);
    temp = sf_floatalloc(n);
    known = sf_boolalloc(n);

    xmov = sf_floatalloc2(n,niter);

    /* specify known locations */
    mask = sf_input("mask");
    sf_floatread(mm,n,mask);
    sf_fileclose(mask);
	
    for (i=0; i < n; i++) {
	known[i] = (bool) (mm[i] != 0.);
	temp[i] = 0.;
    }

    /* Read input data */
    sf_floatread(mm,n,in);

    if (prec) {                          
	sf_mask_init(known);
	sf_solver_prec(sf_mask_lop, sf_cgstep, polydiv_lop, 
		       n, n, n, mm, mm, niter, 0., "end");
    } else {                             
	sf_helicon_init(aa);
	sf_solver (sf_helicon_lop, sf_cgstep, n, n, mm, temp, 
		   niter, "known", known, "x0", mm, "end");
    }

    sf_floatwrite(mm,n,out);

    exit (0);
}
