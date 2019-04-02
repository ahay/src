/* Multidimensional convolution cascade.*/
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

#include "regrid.h"

int main(int argc, char* argv[])
{
    int i, ia, na, ix, nx, dim, nr, ir, n[SF_MAX_DIM], m[SF_MAX_DIM];
    float a0, r, scale, *pp, *qq, *tt;
    sf_filter aa;
    char* lagfile;
    sf_file in, out, filt, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    out = sf_output("out");

    if (!sf_getint("rect",&nr)) nr=0;
    /* smoothing radius */

    dim = sf_filedims (in,n);

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");
    aa = sf_allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
    sf_floatread (aa->flt,na,filt);
    for( ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    if (NULL != (lagfile = sf_getstring("lag")) 
	/*( lag file with filter lags )*/
	|| 
	NULL != (lagfile = sf_histstring(filt,"lag"))) {
	lag = sf_input(lagfile);

	sf_intread(aa->lag,na,lag);
    } else {
	lag = NULL;
	for( ia=0; ia < na; ia++) {
	    aa->lag[ia] = ia+1;
	}
    }

    sf_fileclose(filt);
    
    if (!sf_getints ("n",m,dim) && (NULL == lag ||
				    !sf_histints (lag,"n",m,dim))) {
	for (i=0; i < dim; i++) {
	    m[i] = n[i];
	}
    }
 
    if (NULL != lag) sf_fileclose(lag);

    regrid (dim, m, n, aa);

    nx = 1;
    for( i=0; i < dim; i++) {
	nx *= n[i];
    }
  
    pp = sf_floatalloc (nx);
    qq = sf_floatalloc (nx);
    tt = sf_floatalloc (nx);

    sf_floatread (pp,nx,in);
    sf_helicon_init (aa);

    scale = 1.0f;
    for( ia=0; ia < na; ia++) {
	scale += aa->flt[ia] * aa->flt[ia];
    }
    sf_warning("scale=%g",scale);

    for (ir=1; ir < nr; ir++) {
	r = 2*sinf(SF_PI*ir/nr);
	r = 2/(r*r*scale);
	    
	sf_helicon_lop (false,false,nx,nx,pp,tt);
	sf_helicon_lop (true,false,nx,nx,qq,tt);

	for (ix=0; ix < nx; ix++) {
	    pp[ix] -= r*qq[ix];
	}
    }
 
    sf_floatwrite (pp,nx,out);

    exit (0);
}
