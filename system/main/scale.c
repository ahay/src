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

static float getscale(int ntest, int n1, float *abuf);

int main (int argc, char* argv[])
{
    off_t n1, n2, n[SF_MAX_DIM];
    int axis, ndim, i, i2, *ibuf;
    size_t nsize, nbuf, ntest;
    sf_file in=NULL;
    sf_file out=NULL;
    float *fbuf, *abuf, dscale, pclip;
    sf_complex* cbuf;
    sf_datatype type;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    type = sf_gettype(in);
    nbuf = sf_bufsiz(in);
    if (SF_INT == type) sf_settype(out,SF_FLOAT);

    sf_fileflush(out,in);

    if (!sf_getint("axis",&axis)) axis = 0;
    /* Scale by maximum in the dimensions up to this axis. */
    if (!sf_getfloat("rscale",&dscale)) dscale=0.;
    /* Scale by this factor. */
    if (!sf_getfloat("pclip",&pclip)) pclip=100.;
    /* data clip percentile */

    ndim = sf_largefiledims (in, n);

    n1=1;
    n2=1;

    for (i=0; i < ndim; i++) {
	if (i < axis) n1 *= n[i];
	else          n2 *= n[i];
    }

    abuf = sf_floatalloc(n1);

    ntest = n1*pclip/100. + .5;
    if (ntest < 0)       ntest = 0;
    else if (ntest > n1) ntest = n1;

    if (1 < n1 && 0. == dscale) {
	switch (type) {
	    case SF_FLOAT:
		fbuf = sf_floatalloc(n1);

		for (i2=0; i2 < n2; i2++) {
		    sf_floatread (fbuf,n1,in);
		    for (i=0; i < n1; i++) {
			abuf[i] = fabsf(fbuf[i]);
		    }
		    dscale = getscale(ntest,n1,abuf);		    
		    for (i=0; i < n1; i++) {
			fbuf[i] *= dscale;
		    }
		    sf_floatwrite (fbuf,n1,out);
		}
		break;
	    case SF_COMPLEX:
		cbuf = sf_complexalloc(n1);

		for (i2=0; i2 < n2; i2++) {
		    sf_complexread (cbuf,n1,in);
		    for (i=0; i < n1; i++) {
			abuf[i] = cabsf(cbuf[i]);
		    }
		    dscale = getscale(ntest,n1,abuf);
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

static float getscale(int ntest, int n1, float *abuf)
{
    int i;
    float dscale, f;

    if (ntest == n1) { /* 100% clip - take maximum */
	dscale = 0.;
	for (i=0; i < n1; i++) {
	    f = abuf[i];
	    if (f > dscale) dscale=f;
	}
    } else {
	dscale = sf_quantile(ntest,n1,abuf);
    }

    if (dscale > 0.) dscale=1./dscale;

    return dscale;
}
