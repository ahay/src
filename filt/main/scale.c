/* Scale data.

Takes: < input.rsf > scaled.rsf

To scale by a constant factor, you can also use sfmath of sfheadermath.
*/

#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int axis, n[SF_MAX_DIM], ndim, i, i2, n1, n2, nsize, nbuf = BUFSIZ;
    sf_file in, out;    
    float* fbuf, f, dscale;
    float complex* cbuf;
    sf_datatype type;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");
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
	switch (type = sf_gettype(in)) {
	    case SF_FLOAT:
		fbuf = sf_floatalloc(n1);

		for (i2=0; i2 < n2; i2++) {
		    dscale = 0.;
		    sf_read (fbuf,sizeof(float),n1,in);
		    for (i=0; i < n1; i++) {
			f = fabsf(fbuf[i]);
			if (f > dscale) dscale=f;
		    }
		    if (dscale > 0.) dscale=1./dscale;
		    for (i=0; i < n1; i++) {
			fbuf[i] *= dscale;
		    }
		    sf_write (fbuf,sizeof(float),n1,out);
		}
		break;
	    case SF_COMPLEX:
		cbuf = sf_complexalloc(n1);

		for (i2=0; i2 < n2; i2++) {
		    dscale = 0.;
		    sf_read (cbuf,sizeof(float complex),n1,in);
		    for (i=0; i < n1; i++) {
			f = cabsf(cbuf[i]);
			if (f > dscale) dscale=f;
		    }
		    if (dscale > 0.) dscale=1./dscale;
		    for (i=0; i < n1; i++) {
			cbuf[i] *= dscale;
		    }
		    sf_write (cbuf,sizeof(float complex),n1,out);
		}
		break;
	    default:
		sf_error("Unsupported type %d",type);	
		break;
	}
    } else {
	if (0.==dscale && !sf_getfloat("dscale",&dscale)) dscale=1.;
	/* Scale by this factor (works if rscale=0) */

	nbuf /= sizeof(float);
	fbuf = sf_floatalloc(nbuf);
	nsize = n1*n2;

	switch (type = sf_gettype(in)) {
	    case SF_FLOAT:
		break;
	    case SF_COMPLEX:
		nsize *= 2;
		sf_settype(in,SF_FLOAT);
		sf_settype(out,SF_FLOAT);
		break;
	    default:
		sf_error("Unsupported type %d",type);
	}			
	
	while (nsize > 0) {
	    if (nsize < nbuf) nbuf = nsize;
	    sf_read (fbuf,sizeof(float),nbuf,in);
	    for (i=0; i < nbuf; i++) {
		fbuf[i] *= dscale;
	    }
	    sf_write (fbuf,sizeof(float),nbuf,out);
	    nsize -= nbuf;
	}	
    } 

    sf_close();
    exit (0);
}

/* 	$Id: scale.c,v 1.4 2004/04/06 02:02:54 fomels Exp $	 */
