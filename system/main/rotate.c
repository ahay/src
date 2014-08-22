/* Rotate a portion of one or more axes in the data hypercube. */
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

static void rotate (size_t n1, int dim, 
		    const off_t* n, const int* rot, /*@out@*/ off_t *k);

int main(int argc, char* argv[])
{
    sf_file in, out;
    char *buf, *buf2, key[5];
/* Want them to be arbitrary, neither float nor complex */
/* Just pretend they are character pointers so we multiply offsets ourselves.*/
    int i, dim, dim1, dim2, rot[SF_MAX_DIM], esize;
    int mem; /* for avoiding int to off_t typecast warning */
    off_t n[SF_MAX_DIM], pos=0, pos3=0, memsize, *k1 = NULL, *k2 = NULL;
    size_t n1, i1, i2, i3, n2, n3, size;
    bool verb;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    dim   = sf_largefiledims(in,n);
    esize = sf_esize(in);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* Verbosity flag */

    if (!sf_getint("memsize",&mem))
        mem=sf_memsize();
    /* Max amount of RAM (in Mb) to be used */
    memsize = (off_t) mem * (1<<20); /* convert Mb to bytes */

    dim2=0;
    for (i=0; i < dim; i++) {	
	snprintf(key,5,"rot%d",i+1);
	if (!sf_getint(key,rot+i)) rot[i]=0;
	/*( rot#=(0,0,...) length of #-th axis that is moved to the end )*/
	if (rot[i] < 0) rot[i] += n[i];
	if (rot[i] >= n[i]) 
	    sf_error("rot%d=%d must be smaller than n%d=%d",
		     i,rot[i],i,n[i]);
	if (rot[i]) dim2=i;
    }
    dim2++;
    /* dim2 is the last dimension involved */

    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
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

    if (n1>1) {
	k1 = sf_largeintalloc(n1);
	rotate(n1,dim1,n,rot,k1);
    }

    if (n2>1) {
	if (verb) sf_warning("Going out of core...");

	sf_unpipe(in,n1*n2*n3*esize); /* prepare for random access */
	pos = sf_tell(in);

	k2 = sf_largeintalloc(n2);
	rotate(n2,dim2-dim1,n+dim1,rot+dim1,k2);
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

static void rotate (size_t n1, int dim, 
		    const off_t* n, const int* rot, /*@out@*/ off_t *k)
/* compute map of rotations k */
{
    off_t i, m, ii;
    int j, nj, r1, r2;
    
    for (i=0; i < n1; i++) {
	m=0;
	for (j=0, nj=1; j< dim; nj*=n[j], j++) {
	    ii = (i/nj)%n[j]; /* convert from helix to cartesian */
	    r1 = rot[j]; 
	    r2 = n[j]-r1;
	    if (ii < r2) ii += r1;
	    else         ii -= r2;
	    m += ii*nj;       /* convert from cartesian to helix */
	}
	k[i]=m;
    }
}

/* 	$Id: rotate.c 1729 2006-03-12 10:00:32Z fomels $	 */
