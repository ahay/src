/* Reverse one or more axes in the data hypercube.
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
#include <string.h>

#include <rsf.h>

static void mirror (size_t n1, int dim, 
		    const int* n, const bool* f, /*@out@*/ size_t *k);

int main(int argc, char* argv[])
{
    sf_file in, out;
    char *buf, *opt, copt, key[3], byte;
/* Want them to be arbitrary, neither float nor complex */
/* Just pretend they are character pointers so we multiply offsets ourselves.*/
    int j, i, dim, dim1, dim2;
    int n[SF_MAX_DIM], esize, which, memsize;
    long pos=0;
    size_t n1, i1, i2, i3, n2, n3, size, *k1 = NULL, *k2 = NULL, m;
    unsigned int mask;
    bool f[SF_MAX_DIM], verb;
/* Flags; 0=leave alone, 1=reverse this dimension */
    float o, d;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    dim = sf_filedims(in,n);
    if (!sf_histint(in,"esize",&esize)) esize=4;
    
    if (!sf_getint("which",&which)) which=-1;
    /* Which axis to reverse.
       To reverse a given axis, start with 0,
       add 1 to number to reverse n1 dimension,
       add 2 to number to reverse n2 dimension,
       add 4 to number to reverse n3 dimension, etc.
       Thus, which=7 would reverse the first three dimensions,
       which=5 just n1 and n3, etc.
       which=0 will just pass the input on through unchanged. */

    opt = sf_getstring("opt");
    /* If y, change o and d parameters on the reversed axis. */
    copt = (NULL == opt)? 'y':opt[0];

    /* Figure out which dimension is the slowest */
    if (-1 == which) which = (1 << (dim-1));

    if (!sf_getbool("verb",&verb)) verb=false;
    /* Verbosity flag */
    if (!sf_getint("memsize",&memsize) || 0 >= memsize) 
	memsize = 1 << 20; /* 1 Mb */
   
    if (verb) fprintf(stderr,"%s: Reversing over",sf_getprog());
    for (i=0, mask=1; i < dim; i++, mask <<= 1) {
	f[i] = (0 != (which & mask)) && n[i]>1;
	if (verb && f[i]) fprintf(stderr," n%d",i+1);
    }
    if (verb) fprintf(stderr,".\n");

    dim2=0;
    for (i=0; i < dim; i++) {	
	if (!f[i]) continue;
	dim2=i;
	if ('i' != copt) {     /* update the o's and d's .*/
	    snprintf(key,3,"d%d",i+1);
	    if (!sf_histfloat(in,key,&d)) d=1.;
	    if ('y'==copt) sf_putfloat(out,key,-d);
	    snprintf(key,3,"o%d",i+1); 
	    if (!sf_histfloat(in,key,&o)) o=0.;
	    o += (n[i] - 1) * d;
	    if ('n'==copt) o=-o;
	    sf_putfloat(out,key,o);
	}
    }
    dim2++;

    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    n1=n2=n3=1;
    dim1=0;
    for (i=0; i < dim; i++) {
	if (i < dim2) {
	    size = n1*n[i];
	    if (size*esize < (size_t) memsize) {
		n1=size;
		dim1=i;
	    } else {
		n2*=n[i];
	    }
	} else {
	    n3 *= n[i];
	}
    }
    dim1++;
    
    buf = sf_charalloc(n1*esize);
    
    if (n1>1) {
	k1 = (size_t*) sf_alloc(n1,sizeof(size_t));
	mirror(n1,dim1,n,f,k1);
    }

    if (n2>1) {
	if (verb) sf_warning("Going out of core...");
	sf_unpipe(in,n1*n2*n3*esize);
	k2 = (size_t*) sf_alloc(n2,sizeof(size_t));
	mirror(n2,dim2-dim1,n+dim1,f+dim1,k2);
	pos = sf_tell(in);
    } 

    if (verb) sf_warning("n1=%d, n2=%d, n3=%d",
			 (int) n1,(int) n2,(int) n3);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    if (n2 > 1) sf_seek(in,pos+k2[i2]*n1*esize,SEEK_SET);
	    sf_charread(buf,n1*esize,in);
	    for (i1=0; i1 < n1/2; i1++) {
		m = k1[i1];
		for (j=0; j<esize; j++) {
		    byte = buf[i1*esize+j];
		    buf[i1*esize+j] = buf[m*esize+j];
		    buf[m*esize+j] = byte;
		}
	    }
	    sf_charwrite(buf,n1*esize,out);
	}
    }
    
    if (n2 > 1) sf_close();
    exit (0);
}

static void mirror (size_t n1, int dim, 
		    const int* n, const bool* f, /*@out@*/ size_t *k)
{
    size_t i, m, ii;
    int j, nj;
    
    for (i=0; i < n1/2; i++) {
	m=0;
	for (j=0, nj=1; j< dim; nj*=n[j], j++) {
	    ii = (i/nj)%n[j];
	    if (f[j]) ii = n[j]-1-ii;
	    m += ii*nj;
	}
	k[i]=m;
	k[m]=i;
    }
    if (0 != n1%2) k[n1/2]=n1/2; /* Take care of odd n1 */
}

/* 	$Id$	 */
