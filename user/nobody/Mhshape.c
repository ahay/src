/* Multidimensional shaping using helix transform. */
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

#include "hshape.h"
#include "regrid.h"

int main(int argc, char* argv[])
{
    int i, ia, na, nx, ns, dim, n[SF_MAX_DIM], m[SF_MAX_DIM];
    float a0, *pp, *qq;
    bool adj;
    sf_filter aa;
    char* lagfile;
    sf_file in, out, filt, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    out = sf_output("out");

    dim = sf_filedims (in,n);

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");
    aa = sf_allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
    sf_floatread (aa->flt,na,filt);
    for( ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    if (NULL != (lagfile = sf_getstring("lag")) /* file with filter lags */
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

    if (!sf_getbool ("adj",&adj)) adj=false;
    /* if y, do adjoint operation */
    if (!sf_getint ("ns",&ns)) sf_error("Need ns=");
    /* scaling */

    nx = 1;
    for( i=0; i < dim; i++) {
	nx *= n[i];
    }
  
    pp = sf_floatalloc (nx);
    qq = sf_floatalloc (nx);

    if (adj) {
	sf_floatread (qq,nx,in);
    } else {
	sf_floatread (pp,nx,in);
    }

    hshape_init (nx,ns,aa);
    hshape_lop (adj,false,nx,nx,pp,qq);

    if (adj) {
	sf_floatwrite (pp,nx,out);
    } else {
	sf_floatwrite (qq,nx,out);
    }


    exit (0);
}
