/* Add, multiply, or divide  RSF datasets.

Takes: [< file0.rsf] file1.rsf file2.rsf ...

sfadd operates on integer, float, or complex data, but all the input
and output files must be of the same data type.

An alternative to sfadd is sfmath, which is more versatile, but may be
less efficient.
*/
/*
  Copyright (C) 2006 Colorado School of Mines
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

static void check_compat(int        esize, 
			 size_t     nin, 
			 sf_file*   in, 
			 int        dim, 
			 const int *n);

static void add_float   (bool        collect, 
			 char        cmode, 
			 size_t      nbuf, 
			 float*      buf, 
			 float*      bufi, 
			 float       scale, 
			 float       add);

static void add_int     (bool        collect, 
			 char        cmode, 
			 size_t      nbuf, 
			 int*        buf, 
			 int*        bufi, 
			 float       scale, 
			 float       add);

static void add_complex (bool        collect, 
			 char        cmode, 
			 size_t      nbuf, 
			 sf_complex* buf, 
			 sf_complex* bufi, 
			 float       scale, 
			 float       add);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
int main (int argc, char* argv[])
{
    int      i, dim, n[SF_MAX_DIM], esize;
    size_t   j, nin;
    sf_file *in, out;
    float   *scale, *add;
    bool collect;
    char cmode, *mode;
    sf_datatype type;

    int   mem;
    off_t memsize;
    
    size_t nbuf,nsiz;
    char *bufi,*bufo;

    /* init RSF */
    sf_init (argc, argv);

    /*------------------------------------------------------------*/

    in = (sf_file*) sf_alloc ((size_t) argc,sizeof(sf_file));
    out = sf_output ("out");

    /*------------------------------------------------------------*/

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
    /*------------------------------------------------------------*/
    /* default coefficients and flags */
    scale     = sf_floatalloc (nin);
    add       = sf_floatalloc (nin);      
    for (j = 0; j < nin; j++) {
	scale[j] = 1.;
	add[j]   = 0.;
    }
    (void) sf_getfloats("scale",scale,nin); 
    (void) sf_getfloats("add",add,nin);
    /*------------------------------------------------------------*/
    /* 'a' means add (default), 
       'p' or 'm' means multiply, 
       'd' means divide 
    */
    mode = sf_getstring("mode");
    cmode = (NULL==mode)? 'a':mode[0];
    /*------------------------------------------------------------*/
    /* verify file compatibility */
    dim = sf_filedims(in[0],n); /* input files dimensions */
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }                           /* number of elements in input files */
    esize = (int) sf_esize(in[0]);
    check_compat(esize,nin,in,dim,n);
    /*------------------------------------------------------------*/
    /* prepare output file */
    sf_setformat(out,sf_histstring(in[0],"data_format"));
    sf_fileflush(out,in[0]);
    /*------------------------------------------------------------*/

    if (!sf_getint("memsize",&mem)) mem = 100;
    /* Available memory size (in Mb) */
    memsize = mem * (1 << 20); /* convert Mb to bytes */

    sf_warning("nsiz=%d",nsiz);
    nbuf = memsize/esize/2;
    if(nbuf>nsiz) nbuf=nsiz;
    sf_warning("nbuf=%d",nbuf);

    bufi = (char*) sf_alloc( nbuf, sizeof(char));
    bufo = (char*) sf_alloc( nbuf, sizeof(char));
    
    /*------------------------------------------------------------*/

    type = sf_gettype (out); /* input/output files format */
    
    for (nbuf /= esize; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf=nsiz;

	switch(type) {
	    case SF_FLOAT:
		for (j=0; j < nin; j++) {
		    collect = (bool) (j != 0);
		    sf_floatread((float*) bufi,nbuf,in[j]);	    
		    add_float(collect, cmode, nbuf,
			      (float*) bufo, (float*) bufi, 
			      scale[j], add[j]
			);
		}
		sf_floatwrite((float*) bufo,nbuf,out);
		break;
	    case SF_COMPLEX:
		for (j=0; j < nin; j++) {
		    collect = (bool) (j != 0);
		    sf_complexread((sf_complex*) bufi,nbuf,in[j]);
		    add_complex(collect, cmode, nbuf,
				(sf_complex*) bufo, (sf_complex*) bufi, 
				scale[j], add[j]
			);
		}
		sf_complexwrite((sf_complex*) bufo,nbuf,out);
		break;
	    case SF_INT:
		for (j=0; j < nin; j++) {
		    collect = (bool) (j != 0);
		    sf_intread((int*) bufi,nbuf,in[j]);
		    add_int(collect, cmode, nbuf,
			    (int*) bufo,(int*) bufi,
			    scale[j], add[j]
			);
		}
		sf_intwrite((int*) bufo,nbuf,out);
		break;
	    default:
		sf_error("wrong type");
		break;
	}	
    }
    exit (0);
}

/*------------------------------------------------------------*/
static void add_float (bool   collect, 
		       char   cmode, 
		       size_t nbuf, 
		       float* bufo, 
		       float* bufi, 
		       float  scale, 
		       float  add)
{
    size_t jbuf;
    float f;

    if (collect) {
	switch (cmode) {
	    case 'p':
	    case 'm':
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    bufo[jbuf] *= scale*(bufi[jbuf] + add);
		}
		break;
	    case 'd':
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf,f) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    f = scale*(bufi[jbuf] + add);
		    if (f != 0.) bufo[jbuf] /= f;
		}
		break;
	    default:
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    bufo[jbuf] +=scale*(bufi[jbuf] + add);
		}
		break;
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf) shared(scale,add,bufi) 
#endif
	for(jbuf=0;jbuf<nbuf;jbuf++){
	    bufo[jbuf] = scale*(bufi[jbuf] + add);
	}
    }
}

/*------------------------------------------------------------*/
static void add_int (bool   collect, 
		     char   cmode, 
		     size_t nbuf, 
		     int*   bufo, 
		     int*   bufi, 
		     float  scale, 
		     float  add)
{
    size_t jbuf;
    float f;

    if (collect) {
	switch (cmode) {
	    case 'p':
	    case 'm':
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    bufo[jbuf] *= scale*( (float) bufi[jbuf] + add);
		}
		break;
	    case 'd':
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf,f) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    f = scale*( (float) bufi[jbuf] + add);
		    if (f != 0.) bufo[jbuf] /= f;
		}
		break;
	    default:
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    bufo[jbuf] +=scale*( (float) bufi[jbuf] + add);
		}
		break;
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) private(jbuf) shared(scale,add,bufi) 
#endif
	for(jbuf=0;jbuf<nbuf;jbuf++){
	    bufo[jbuf] = scale*( (float) bufi[jbuf] + add);
	}
    }
}

/*------------------------------------------------------------*/
static void add_complex (bool        collect, 
			 char        cmode, 
			 size_t      nbuf, 
			 sf_complex* bufo, 
			 sf_complex* bufi, 
			 float       scale, 
			 float       add)
{
    size_t j;
    sf_complex c;
    
    for (j=0; j < nbuf; j++) {
	c = bufi[j];
#ifdef SF_HAS_COMPLEX_H
	c += add;
#else
	c.r += add;
#endif
	if (1. != scale) {
#ifdef SF_HAS_COMPLEX_H
	    c *= scale;
#else
	    c = sf_crmul(c,scale);
#endif
	}
	if (collect) {
	    switch (cmode) {
		case 'p':
		case 'm':
#ifdef SF_HAS_COMPLEX_H
		    bufo[j] *= c;
#else
		    bufo[j] = sf_cmul(bufo[j],c);
#endif
		    break;
		case 'd':
		    if (cabsf(c) != 0.) {
#ifdef SF_HAS_COMPLEX_H			
			bufo[j] /= c;
#else
			bufo[j] = sf_cdiv(bufo[j],c);
#endif	
		    }
		    break;
		default:
#ifdef SF_HAS_COMPLEX_H	
		    bufo[j] += c;
#else
		    bufo[j] = sf_cadd(bufo[j],c);
#endif
		    break;
	    }
	} else {
	    bufo[j] = c;
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

