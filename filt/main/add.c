/* Add, multiply, or divide  RSF datasets.

Takes: [< file0.rsf] file1.rsf file2.rsf ...

The various operations, if selected, occur in the following order:

(1) Take absolute value, abs=
(2) Add a scalar, add=
(3) Take the natural logarithm, log=
(4) Take the square root, sqrt=
(5) Multiply by a scalar, scale=
(6) Compute the base-e exponential, exp=
(7) Add, multiply, or divide the data sets, mode=

sfadd operates on integer, float, or complex data, but all the input
and output files must be of the same data type.

An alternative to sfadd is sfmath, which is more versatile, but may be
less efficient.
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
#include <string.h>

#include <unistd.h>

#include <rsf.h>

static void check_compat(int        esize, 
			 size_t     nin, 
			 sf_file*   in, 
			 int        dim, 
			 const int *n);
static void add_float (bool   collect, 
		       size_t nbuf, 
		       float* buf, 
		       float* bufi, 
		       char   cmode, 
		       float  scale, 
		       float  add, 
		       bool   abs_flag, 
		       bool   log_flag, 
		       bool   sqrt_flag, 
		       bool   exp_flag);
static void add_int (bool     collect, 
		     size_t   nbuf, 
		     int*     buf, 
		     int*     bufi, 
		     char     cmode, 
		     float    scale, 
		     float    add, 
		     bool     abs_flag, 
		     bool     log_flag, 
		     bool     sqrt_flag, 
		     bool     exp_flag);
static void add_complex (bool        collect, 
			 size_t      nbuf, 
			 sf_complex* buf, 
			 sf_complex* bufi, 
			 char        cmode, 
			 float       scale, 
			 float       add, 
			 bool        abs_flag, 
			 bool        log_flag, 
			 bool        sqrt_flag, 
			 bool        exp_flag);

int main (int argc, char* argv[])
{
    int i, dim, n[SF_MAX_DIM], esize;
    size_t j, nin, nbuf, nsiz;
    sf_file *in, out;
    float *scale, *add;
    bool *sqrt_flag, *abs_flag, *log_flag, *exp_flag, collect;
    char cmode, *mode, *buf, *bufi;
    sf_datatype type;

    /* init RSF */
    sf_init (argc, argv);

    in = (sf_file*) sf_alloc ((size_t) argc,sizeof(sf_file));
    out = sf_output ("out");
    
    /* find number of input files */
    if (0 != isatty(fileno(stdin))) { /* no input file in stdin */
	nin=0;
    } else {
	in[0] = sf_input("in");
	nin=1;
    }

    for (i=1; i< argc; i++) { /* collect inputs */
	if (NULL != strchr(argv[i],'=')) continue; /* not a file */
	in[nin] = sf_input(argv[i]);
	nin++;
    }
    if (0==nin) sf_error ("no input");
    /* nin = no of input files*/

    nbuf = sf_bufsiz(in[0]);
    buf = sf_charalloc(nbuf);
    bufi = sf_charalloc(nbuf);
    
    out = sf_output ("out");
    
    scale = sf_floatalloc (nin);
    add   = sf_floatalloc (nin);  

    /* default coefficients and flags */
    scale     = sf_floatalloc (nin);
    add       = sf_floatalloc (nin);  
    sqrt_flag = sf_boolalloc  (nin);
    abs_flag  = sf_boolalloc  (nin);
    log_flag  = sf_boolalloc  (nin);
    exp_flag  = sf_boolalloc  (nin);
    
    for (j = 0; j < nin; j++) {
	scale[j] = 1.;
	add[j]   = 0.;
	sqrt_flag[j] = abs_flag[j] = log_flag[j] = exp_flag[j] = false;
    }

    /* Scalar values to multiply each dataset with */
    (void) sf_getfloats("scale",scale,nin); 

    /* Scalar values to add to each dataset */
    (void) sf_getfloats("add",add,nin);

    /* If true (1) take square root */
    (void) sf_getbools("sqrt",sqrt_flag,nin);
    /* If true (1) take absolute value */
    (void) sf_getbools("abs",abs_flag,nin);
    /* If true (1) take logarithm */
    (void) sf_getbools("log",log_flag,nin);
    /* If true (1) compute exponential */
    (void) sf_getbools("exp",exp_flag,nin);

    /* 'a' means add (default), 
       'p' or 'm' means multiply, 
       'd' means divide 
    */
    mode = sf_getstring("mode");
    cmode = (NULL==mode)? 'a':mode[0];

    /* verify file compatibility */
    dim = sf_filedims(in[0],n); /* input files dimensions */
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }                           /* number of elements in input files */
    esize = (int) sf_esize(in[0]);
    check_compat(esize,nin,in,dim,n);
    /* end verify file compatibility */

    sf_setformat(out,sf_histstring(in[0],"data_format"));
    sf_fileflush(out,in[0]);

    type = sf_gettype (out); /* input/output files format */
    
    for (nbuf /= esize; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf=nsiz;

	for (j=0; j < nin; j++) {
	    collect = (bool) (j != 0);
	    switch(type) {
		case SF_FLOAT:
		    sf_floatread((float*) bufi,
				 nbuf,
				 in[j]);	    
		    add_float(collect, 
			      nbuf,
			      (float*) buf,
			      (float*) bufi, 
			      cmode, 
			      scale[j], 
			      add[j], 
			      abs_flag[j], 
			      log_flag[j], 
			      sqrt_flag[j], 
			      exp_flag[j]);
		    break;
		case SF_COMPLEX:		    
		    sf_complexread((sf_complex*) bufi,
				   nbuf,
				   in[j]);
		    add_complex(collect, 
				nbuf,
				(sf_complex*) buf,
				(sf_complex*) bufi, 
				cmode, 
				scale[j], 
				add[j], 
				abs_flag[j], 
				log_flag[j], 
				sqrt_flag[j], 
				exp_flag[j]);
		    break;
		case SF_INT:
		    sf_intread((int*) bufi,
			       nbuf,
			       in[j]);
		    add_int(collect, 
			    nbuf,
			    (int*) buf,
			    (int*) bufi, 
			    cmode, 
			    scale[j], 
			    add[j], 
			    abs_flag[j], 
			    log_flag[j], 
			    sqrt_flag[j], 
			    exp_flag[j]);
		    break;
		default:
		    sf_error("wrong type");
		    break;
	    }
	}

	switch(type) {
	    case SF_FLOAT:
		sf_floatwrite((float*) buf,nbuf,out);
		break;
	    case SF_COMPLEX:
		sf_complexwrite((sf_complex*) buf,nbuf,out);
		break;
	    case SF_INT:
		sf_intwrite((int*) buf,nbuf,out);
		break;
	    default:
		sf_error("wrong type");
		break;
	}	
	
    }
    
    exit (0);
}

static void add_float (bool   collect, 
		       size_t nbuf, 
		       float* buf, 
		       float* bufi, 
		       char   cmode, 
		       float  scale, 
		       float  add, 
		       bool   abs_flag, 
		       bool   log_flag, 
		       bool   sqrt_flag, 
		       bool   exp_flag)
{
    size_t j;
    float f;

    for (j=0; j < nbuf; j++) {
	f = bufi[j];
	if (abs_flag)    f = fabsf(f);
	f += add;
	if (log_flag)    f = logf(f);
	if (sqrt_flag)   f = sqrtf(f);
	if (1. != scale) f *= scale;
	if (exp_flag)    f = expf(f);
	if (collect) {
	    switch (cmode) {
		case 'p':
		case 'm':
		    buf[j] *= f;
		    break;
		case 'd':
		    if (f != 0.) buf[j] /= f;
		    break;
		default:
		    buf[j] += f;
		    break;
	    }
	} else {
	    buf[j] = f;
	}
    }
}

static void add_int (bool   collect, 
		     size_t nbuf, 
		     int*   buf, 
		     int*   bufi, 
		     char   cmode, 
		     float  scale, 
		     float  add, 
		     bool   abs_flag, 
		     bool   log_flag, 
		     bool   sqrt_flag, 
		     bool   exp_flag)
{
    size_t j;
    float f;

    for (j=0; j < nbuf; j++) {
	f = (float) bufi[j];
	if (abs_flag)    f = fabsf(f);
	f += add;
	if (log_flag)    f = logf(f);
	if (sqrt_flag)   f = sqrtf(f);
	if (1. != scale) f *= scale;
	if (exp_flag)    f = expf(f);
	if (collect) {
	    switch (cmode) {
		case 'p':
		case 'm':
		    buf[j] *= f;
		    break;
		case 'd':
		    if (f != 0.) buf[j] /= f;
		    break;
		default:
		    buf[j] += f;
		    break;
	    }
	} else {
	    buf[j] = f;
	}
    }
}

static void add_complex (bool        collect, 
			 size_t      nbuf, 
			 sf_complex* buf, 
			 sf_complex* bufi, 
			 char        cmode, 
			 float       scale, 
			 float       add, 
			 bool        abs_flag, 
			 bool        log_flag, 
			 bool        sqrt_flag, 
			 bool        exp_flag)
{
    size_t j;
    sf_complex c;
    
    for (j=0; j < nbuf; j++) {
	c = bufi[j];
	if (abs_flag)  {
#ifdef SF_HAS_COMPLEX_H
	    c = cabsf(c);
#else
	    c.r = cabsf(c);
	    c.i = 0.;
#endif
	}
#ifdef SF_HAS_COMPLEX_H
	c += add;
#else
	c.r += add;
#endif
	if (log_flag)    c = clogf(c);
	if (sqrt_flag)   c = csqrtf(c);
	if (1. != scale) {
#ifdef SF_HAS_COMPLEX_H
	    c *= scale;
#else
	    c = sf_crmul(c,scale);
#endif
	}
	if (exp_flag)    c = cexpf(c);
	if (collect) {
	    switch (cmode) {
		case 'p':
		case 'm':
#ifdef SF_HAS_COMPLEX_H
		    buf[j] *= c;
#else
		    buf[j] = sf_cmul(buf[j],c);
#endif
		    break;
		case 'd':
		    if (cabsf(c) != 0.) {
#ifdef SF_HAS_COMPLEX_H			
			buf[j] /= c;
#else
			buf[j] = sf_cdiv(buf[j],c);
#endif	
		    }
		    break;
		default:
#ifdef SF_HAS_COMPLEX_H	
		    buf[j] += c;
#else
		    buf[j] = sf_cadd(buf[j],c);
#endif
		    break;
	    }
	} else {
	    buf[j] = c;
	}
    }
}

/*------------------------------------------------------------*/
static void check_compat (int        esize, 
			  size_t     nin, 
			  sf_file*   in, 
			  int        dim, 
			  const int* n) 
{
    int ni, id;
    size_t i;
    float d, di, o, oi;
    char key[3];
    const float tol=1.e-5;
    
    for (i=1; i < nin; i++) {
	if ((int) sf_esize(in[i]) != esize) 
	    sf_error ("esize mismatch: need %d",esize);
	for (id=1; id <= dim; id++) {
	    (void) snprintf(key,3,"n%d",id);
	    if (!sf_histint(in[i],key,&ni) || ni != n[id-1])
		sf_error("%s mismatch: need %d",key,n[id-1]);
	    (void) snprintf(key,3,"d%d",id);
	    if (sf_histfloat(in[0],key,&d)) {
		if (!sf_histfloat(in[i],key,&di) || 
		    (fabsf(di-d) > tol*fabsf(d)))
		    sf_warning("%s mismatch: need %g",key,d);
	    } else {
		d = 1.;
	    }
	    (void) snprintf(key,3,"o%d",id);
	    if (sf_histfloat(in[0],key,&o) && 
		(!sf_histfloat(in[i],key,&oi) || (fabsf(oi-o) > tol*fabsf(d))))
		sf_warning("%s mismatch: need %g",key,o);
	}
    }
}

