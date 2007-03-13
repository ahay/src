/* Scale data.

To scale by a constant factor, you can also use sfmath.
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

#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int axis, n[SF_MAX_DIM], ndim, i, i2, n1, n2, nsize, nbuf = BUFSIZ, *ibuf;
    sf_file in, out;    
    float* fbuf, f, dscale;
    sf_complex* cbuf;
    sf_datatype type;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");
  
    type = sf_gettype(in);
    if (SF_INT == type) sf_settype(out,SF_FLOAT);

    sf_fileflush(out,in);

    if (!sf_getint("axis",&axis)) axis = 0;
    /* Scale by maximum in the dimensions up to this axis. */
    if (!sf_getfloat("rscale",&dscale)) dscale=0.;
    /* Scale by this factor. */

    ndim = sf_filedims (in, n);

    n1=1;
    n2=1;

    for (i=0; i < ndim; i++) {
	if (i < axis) n1 *= n[i];
	else          n2 *= n[i];
    }


    if (1 < n1 && 0. == dscale) {
	switch (type) {
	    case SF_FLOAT:
		fbuf = sf_floatalloc(n1);

		for (i2=0; i2 < n2; i2++) {
		    dscale = 0.;
		    sf_floatread (fbuf,n1,in);
		    for (i=0; i < n1; i++) {
			f = fabsf(fbuf[i]);
			if (f > dscale) dscale=f;
		    }
		    if (dscale > 0.) dscale=1./dscale;
		    for (i=0; i < n1; i++) {
			fbuf[i] *= dscale;
		    }
		    sf_floatwrite (fbuf,n1,out);
		}
		break;
	    case SF_COMPLEX:
		cbuf = sf_complexalloc(n1);

		for (i2=0; i2 < n2; i2++) {
		    dscale = 0.;
		    sf_complexread (cbuf,n1,in);
		    for (i=0; i < n1; i++) {
			f = cabsf(cbuf[i]);
			if (f > dscale) dscale=f;
		    }
		    if (dscale > 0.) dscale=1./dscale;
		    for (i=0; i < n1; i++) {
#ifdef SF_HAS_COMPLEX_H
			cbuf[i] *= dscale;
#else
			cbuf[i] = sf_crmul(cbuf[i],dscale);
#endif
		    }
		    sf_complexwrite (cbuf,n1,out);
		}
		break;
	    default:
		sf_error("Unsupported type %d",type);	
		break;
	}
    } else {
	if (0.==dscale && !sf_getfloat("dscale",&dscale)) dscale=1.;
	/* Scale by this factor (works if rscale=0) */

	nsize = n1*n2;

	switch (type) {
	    case SF_COMPLEX:
		nsize *= 2;
		sf_settype(in,SF_FLOAT);
		sf_settype(out,SF_FLOAT);
	    case SF_FLOAT:
		nbuf /= sizeof(float);
		fbuf = sf_floatalloc(nbuf);
	
		while (nsize > 0) {
		    if (nsize < nbuf) nbuf = nsize;
		    sf_floatread (fbuf,nbuf,in);
		    for (i=0; i < nbuf; i++) {
			fbuf[i] *= dscale;
		    }
		    sf_floatwrite (fbuf,nbuf,out);
		    nsize -= nbuf;
		}	
		break;
	    case SF_INT:
		nbuf /= sizeof(int);
		ibuf = sf_intalloc(nbuf);
		fbuf = sf_floatalloc(nbuf);
		
		while (nsize > 0) {
		    if (nsize < nbuf) nbuf = nsize;
		    sf_intread (ibuf,nbuf,in);
		    for (i=0; i < nbuf; i++) {
			fbuf[i] = (double) ibuf[i]*dscale;
		    }
		    sf_floatwrite (fbuf,nbuf,out);
		    nsize -= nbuf;
		}	
		break;
	    default:
		sf_error("Unsupported type %d",type);

	} 
    }

    exit (0);
}

/* 	$Id$	 */
