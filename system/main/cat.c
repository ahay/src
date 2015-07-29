/* Concatenate datasets. 

Takes: [<file0.rsf] file1.rsf file2.rsf ... 

sfmerge inserts additional space between merged data.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <unistd.h>

#include <rsf.h>

static int *order, *sort;

static int sort_order(const void *a, const void *b) 
{
    int ia, ib;

    ia = sort[(int*) a - order];
    ib = sort[(int*) b - order];

    return (ia < ib)? -1: (ia > ib)? 1: 0;
}

static void check_compat (int esize, int nin, int nopen, sf_file *ins, 
			  const char **filename, int axis, 
			  int dim, const off_t *n, /*@out@*/ int *naxis);

int main (int argc, char* argv[])
{
    int i, j, k, axis, *naxis, nin, open_max, nopen;
    int dim, dim1, esize, nspace;
    float f;
    off_t ni, nbuf, n1, n2, i2, n[SF_MAX_DIM], *tell; 
    sf_file *ins, in, out;
    char *prog, key[3], buf[BUFSIZ];
    const char **filename;
    bool space;
    
    sf_init(argc,argv);
    
    filename = (const char**) sf_alloc ((size_t) argc,sizeof(char*));

    if (!sf_stdin()) { /* no input file in stdin */
	nin=0;
    } else {
	filename[0] = "in";
	nin=1;
    }

    for (i=1; i< argc; i++) { /* collect inputs */
	if (NULL != strchr(argv[i],'=')) continue; /* not a file */
	filename[nin] = argv[i];
	nin++;
    }
    if (0==nin) sf_error ("no input");

    order = sf_intalloc(nin);
    for (j=0; j < nin; j++) order[j]=j;

    sort = sf_intalloc(nin);
    if (sf_getints("order",sort,nin)) {
	/* concatenation order */
	qsort(order,nin,sizeof(int),sort_order);
    }
    free(sort);

    open_max = sysconf(_SC_OPEN_MAX);
    /* system limit for the number of simultaneously open files */

    nopen = (open_max > 0)? open_max/2-10:10;

    if (nin > nopen) {
	tell = (off_t*) sf_alloc((size_t) nin-nopen,sizeof(off_t));
	for (i=0; i< nin-nopen; i++) {
	    tell[i] = 0;
	}
    } else {
	tell = NULL;
	nopen = nin;
    }

    ins = (sf_file*) sf_alloc((size_t) nopen,sizeof(sf_file));

    for (i=0; i < nopen; i++) {
	ins[i] = sf_input(filename[i]);
    }
    out = sf_output ("out");

    if (!sf_getbool("space",&space)) {
	/* Insert additional space.
	   y is default for sfmerge, n is default for sfcat */
	prog = sf_getprog();
	if (NULL != strstr (prog, "merge")) {
	    space = true;
	} else if (NULL != strstr (prog, "cat") || 
		   NULL != strstr (prog, "rcat")) {
	    space = false;
	} else {
	    sf_warning("%s is neither merge nor cat,"
		       " assume merge",prog);
	    space = true;
	}
    }

    dim = sf_largefiledims(ins[0],n);
    if (!sf_getint("axis",&axis)) axis=3;
    /* Axis being merged */
    if (0 == axis) {
       axis = SF_MAX_DIM;
       /* Find last axis */
       while (n[axis - 1] <= 1 && axis >= 1)
           axis--;
    }
    if (1 > axis) sf_error("axis=%d < 1",axis);

    dim1 = dim;
    if (axis > dim) {
	while (dim < axis) {
	    n[dim++] = 1;
	}
    }

    n1=1;
    n2=1;
    for (i=1; i <= dim; i++) {
	if      (i < axis) n1 *= n[i-1];
	else if (i > axis) n2 *= n[i-1];
    }

    naxis = sf_intalloc(nin);

    if (!sf_histint(ins[0],"esize",&esize)) {
	esize=4;
    } else if (0>=esize) {
	sf_error("cannot handle esize=%d",esize);
    }
    esize = sf_esize(ins[0]);
    check_compat(esize,nin,nopen,ins,filename,axis,dim1,n,naxis);

    /* figure out the length of extended axis */
    ni = 0;
    for (j=0; j < nin; j++) {
	ni += naxis[j];
    }

    if (space) {
	if (!sf_getint("nspace",&nspace)) 
	    nspace = (int) (ni/(20*nin) + 1);
	/* if space=y, number of traces to insert */ 
	ni += nspace*(nin-1);
    } 

    (void) snprintf(key,3,"n%d",axis);
    sf_putint(out,key,(int) ni);
    
    if (sf_getfloat("o",&f)) {
	/* axis origin */
	(void) snprintf(key,3,"o%d",axis);
	sf_putfloat(out,key,f);
    }

    if (sf_getfloat("d",&f)) {
	/* axis sampling */
	(void) snprintf(key,3,"d%d",axis);	
	sf_putfloat(out,key,f);
    }

    sf_setformat(out,sf_histstring(ins[0],"data_format"));
    sf_fileflush(out,ins[0]);
    sf_setform(out,SF_NATIVE);

    for (i=0; i < nopen; i++) {
	sf_setform(ins[i],SF_NATIVE);
    }

    for (i2=0; i2 < n2; i2++) {
	for (j=0; j < nin; j++) {
	    k = order[j];
	    
	    if (k < nopen) {
		in = ins[k];
	    } else {
		in = sf_input(filename[k]);
		sf_setform(in,SF_NATIVE);
		sf_seek(in,tell[k-nopen],SEEK_SET);
	    }    

	    for (ni = n1*naxis[k]*esize; ni > 0; ni -= nbuf) {
		nbuf = (BUFSIZ < ni)? BUFSIZ: ni;
		sf_charread (buf,nbuf,in);
		sf_charwrite (buf,nbuf,out);
	    }
	    
	    if (k >= nopen) {
		tell[k-nopen] = sf_tell(in);
		sf_fileclose(in);
	    }

	    if (!space || j == nin-1) continue;
	    /* Add spaces */
	    memset(buf,0,BUFSIZ);
	    for (ni = n1*nspace*esize; ni > 0; ni -= nbuf) {
		nbuf = (BUFSIZ < ni)? BUFSIZ: ni;
		sf_charwrite (buf,nbuf,out);
	    }
	}
    }
    
    exit(0);
}

static void check_compat (int esize, int nin, int nopen, sf_file *ins, 
			  const char **filename, int axis, int dim, 
			  const off_t *n, /*@out@*/ int *naxis) 
/* check if the file dimensions are compatible */
{
    int i, ni, id;
    float o, d, di, oi;
    char key[3];
    sf_file in;
    const float tol=1.e-3;
    
    naxis[0] = n[axis-1];
    for (i=1; i < nin; i++) {
	in = (i >= nopen)? sf_input(filename[i]): ins[i];

	if (!sf_histint(in,"esize",&ni) || ni != esize)
	    sf_error ("esize mismatch: need %d",esize);
	for (id=1; id <= dim; id++) {
	    (void) snprintf(key,3,"n%d",id);
	    if (!sf_histint(in,key,&ni) || (id != axis && ni != n[id-1]))
#if defined(__cplusplus) || defined(c_plusplus)
		sf_error("%s mismatch: need %ld",key,(long int) n[id-1]);
#else
	        sf_error("%s mismatch: need %lld",key,(long long int) n[id-1]);
#endif
	    if (id == axis) naxis[i] = ni;
	    (void) snprintf(key,3,"d%d",id);
	    if (sf_histfloat(ins[0],key,&d)) {
		if (!sf_histfloat(in,key,&di) || 
		    (id != axis && fabsf(di-d) > tol*fabsf(d)))
		    sf_warning("%s mismatch: need %g",key,d);
	    } else {
		d = 1.;
	    }
	    (void) snprintf(key,3,"o%d",id);
	    if (sf_histfloat(ins[0],key,&o) && 
		(!sf_histfloat(in,key,&oi) || 
		 (id != axis && fabsf(oi-o) > tol*fabsf(d))))
		sf_warning("%s mismatch: need %g",key,o);
	}
	if (axis > dim) {
	    (void) snprintf(key,3,"n%d",axis);
	    if (!sf_histint(in,key,naxis+i)) naxis[i]=1;
	}
	
	if (i >= nopen) sf_fileclose(in);
    } /* i */
}

/* 	$Id: cat.c 13759 2015-01-09 03:51:11Z sfomel $	 */

