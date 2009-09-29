/* Reverse one or more axes in the data hypercube. */
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

static void mirror (size_t n1, 
		    int dim, 
		    const off_t* n, 
		    const bool* f, 
		    /*@out@*/ off_t *k);

int main(int argc, char* argv[])
{
    sf_file in, out;
    char *buf, *buf2, *opt, copt, key[3];
/* Want them to be arbitrary, neither float nor complex */
/* Just pretend they are character pointers so we multiply offsets ourselves.*/
    int i, dim, dim1, dim2, esize, which;
    int mem; /* for avoiding int to off_t typecast warning */
    off_t n[SF_MAX_DIM], pos=0, pos3=0, memsize, size, *k1 = NULL, *k2 = NULL;
    size_t n1, n2, n3, i1, i2, i3;
    unsigned int mask;
    bool f[SF_MAX_DIM], verb;
/* Flags; 0=leave alone, 1=reverse this dimension */
    float o, d;

    sf_init(argc,argv);

    in  = sf_input(  "in");
    out = sf_output("out");

    dim   = sf_largefiledims(in,n);
    esize = sf_esize(in);
    
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
    /* If y, change o and d parameters on the reversed axis;
       if i, don't change o and d */
    copt = (NULL == opt)? 'y':opt[0];

    /* Figure out which dimension is the slowest */
    if (-1 == which) which = (1 << (dim-1));

    if (!sf_getbool("verb",&verb)) verb=false;
    /* Verbosity flag */

    if (!sf_getint("memsize",&mem))
        mem=sf_memsize();
    /* Max amount of RAM (in Mb) to be used */
    memsize = mem * (1<<20); /* convert Mb to bytes */

    if (verb) fprintf(stderr,"%s: Reversing over",sf_getprog());
    for (i=0, mask=1; i < dim; i++, mask <<= 1) {
	f[i] = (bool) ((0 != (which & mask)) && n[i]>1);
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
    /* dim2 is the last dimension involved */

    sf_fileflush(out,in);
    sf_setform( in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    n1=n2=n3=1;
    dim1=0;
    for (i=0; i < dim; i++) {
	if (i < dim2) {
	    size = n1*n[i];
	    if (1==n2 && size*esize < memsize) {
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
    /* dim1 is the last dimension that fits in memory,
       n1 is the size up to this dimension
       n2 is the size of other dimensions involved 
       n3 is the size of all other dimensions so that
       n1*n2*n3 is the total data size
    */
    
    buf  = sf_charalloc(n1*esize);
    buf2 = sf_charalloc(n1*esize);

    if (n1>=1) {
	k1 = sf_largeintalloc(n1);
	mirror(n1,dim1,n,f,k1);
    }

    if (n2>1) {
	if (verb) sf_warning("Going out of core...");

	sf_unpipe(in,(off_t) n1*n2*n3*esize); /* prepare for random access */
	pos = sf_tell(in);

	k2 = sf_largeintalloc(n2);
	mirror(n2,dim2-dim1,n+dim1,f+dim1,k2);
    } 

    if (verb) sf_warning("n1=%d, n2=%d, n3=%d",
			 (int) n1,(int) n2,(int) n3);

    /* k1 is a table for in-core     reversal 
       k2 is a table for out-of-core reversal */

    for (i3=0; i3 < n3; i3++) {
	if (n2 > 1) pos3 = pos+(off_t) i3*n2*n1*esize;
	for (i2=0; i2 < n2; i2++) {
	    if (n2 > 1) /* if out of core */
		sf_seek(in,pos3+k2[i2]*n1*esize,SEEK_SET);

	    sf_charread(buf,n1*esize,in);
	    for (i1=0; i1 < n1; i1++) {
		memcpy(buf2+k1[i1]*esize,buf+i1*esize,esize);
	    }
	    sf_charwrite(buf2,n1*esize,out);
	}
    }
    

    exit (0);
}

/*------------------------------------------------------------*/

static void mirror (size_t n1, int dim, 
		    const off_t* n, const bool* f, /*@out@*/ off_t *k)
/* compute map of reversals k */
{
    off_t i, m, ii, nj;
    int j;
    
    for (i=0; i < n1; i++) {
	m=0;
	for (j=0, nj=1; j< dim; nj*=n[j], j++) {
	    ii = (i/nj)%n[j]; /* convert from helix to cartesian */
	    if (f[j]) ii = n[j]-1-ii; /* reversal */
	    m += ii*nj;       /* convert from cartesian to helix */
	}
	k[i]=m;
    }
}
