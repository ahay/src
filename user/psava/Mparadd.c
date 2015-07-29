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
			 off_t       nbuf, 
			 float*      buf, 
			 float*      bufi, 
			 float       scale, 
			 float       add,
                         int ompchunk);

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
/*    sf_datatype type; */

    int   mem;
    
    off_t nbuf,nsiz;
    float *bufi,*bufo;
    
    int ompchunk;  /* OpenMP data chunk size */

    /* init RSF */
    sf_init (argc, argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;

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
    /* convert Mb to bytes */
    nbuf = mem/esize*1024*1024;
    if(nsiz>0 && nbuf>nsiz) nbuf=nsiz;
    sf_warning("nbuf=%ld",nbuf);
    
    bufi = sf_floatalloc( nbuf);
    bufo = sf_floatalloc( nbuf);
    
    /*------------------------------------------------------------*/
    /*  type = sf_gettype (out); input/output files format */
    
    for (; nsiz > 0; nsiz -= nbuf) {
	sf_warning("nsiz=%ld nbuf=%ld",nsiz,nbuf);
	if (nbuf > nsiz) nbuf=nsiz;

	for (j=0; j < nin; j++) {
	    collect = (bool) (j != 0);
	    
	    sf_floatread(bufi,nbuf,in[j]);	    
	    
	    add_float(collect, cmode, nbuf,
		      bufo, bufi, 
		      scale[j], add[j],
		      ompchunk
		);
	}
	sf_floatwrite(bufo,nbuf,out);
    }

    free(bufi);
    free(bufo);

    exit (0);
}

/*------------------------------------------------------------*/
static void add_float (bool   collect, 
		       char   cmode, 
		       off_t  nbuf, 
		       float* bufo, 
		       float* bufi, 
		       float  scale, 
		       float  add,
                       int    ompchunk)
{
    int jbuf;
    float f;

    if (collect) {
	switch (cmode) {
	    case 'p':
	    case 'm':
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(jbuf) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    bufo[jbuf] *= scale*(bufi[jbuf] + add);
		}
		break;
	    case 'd':
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(jbuf,f) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    f = scale*(bufi[jbuf] + add);
		    if (f != 0.) bufo[jbuf] /= f;
		}
		break;
	    default:
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(jbuf) shared(scale,add,bufi) 
#endif
		for(jbuf=0;jbuf<nbuf;jbuf++){
		    bufo[jbuf] +=scale*(bufi[jbuf] + add);
		}
		break;
	}
    } else {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(jbuf) shared(scale,add,bufi) 
#endif
	for(jbuf=0;jbuf<nbuf;jbuf++){
	    bufo[jbuf] = scale*(bufi[jbuf] + add);
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

