/* Transpose two axes in a dataset. 

If you get a "Cannot allocate memory" error, give the program a
memsize=1 command-line parameter to force out-of-core operation.
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

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <limits.h>
#include <unistd.h>

#include <stdio.h>
#include <string.h>

#include <rsf.h>

static off_t make_map (int dim1, int dim2, const off_t* n, off_t i2);

int main(int argc, char* argv[])
{
    int i, dim, n3;
    int dim1, dim2;
    int mem; /* for avoiding int to off_t typecast warning */
    off_t n[SF_MAX_DIM], pos, memsize, n1, n2, nsiz, i2, i3, *map, nbuf, n12;
    char key1[7], key2[7], *val, **dat1, **dat2, *buf, *mapf;
    sf_file in, out;
    FILE *mapfile;
    float f;

    sf_init (argc,argv);
    in  = sf_input  ( "in");
    out = sf_output ("out");

    if (!sf_getint("memsize",&mem))
        mem=sf_memsize();
    /* Max amount of RAM (in Mb) to be used */
    memsize = (off_t) mem * (1<<20); /* convert Mb to bytes */

    dim = sf_largefiledims(in,n);

    if (!sf_getint("plane",&dim1)) {
	/* Two-digit number with axes to transpose. The default is 12 */
	dim1=1;
	dim2=2;
    } else {
	dim2 = dim1%10;
	dim1 /= 10;
	
	if (dim2 < 1 || dim2 > SF_MAX_DIM ||
	    dim1 < 1 || dim1 > SF_MAX_DIM || 
	    dim1 == dim2) 
	    sf_error("plane=%d%d is out of range",dim1,dim2);
	
	if (dim1 > dim2) { /* swap so that dim1 is smaller */
	    i = dim1;
	    dim1 = dim2;
	    dim2 = i;
	}
    }
	
    if (dim2 > dim) { /* increase dimensionality if needed */
	for (i=dim; i < dim2; i++) {
	    n[i]=1;
	}
	dim = dim2;
    }

    /* swap axis */
    snprintf(key1,3,"n%d",dim1);
    snprintf(key2,3,"n%d",dim2%10u);
    sf_putint(out,key1,n[dim2-1]);
    sf_putint(out,key2,n[dim1-1]);
 
    snprintf(key1,3,"d%d",dim1);
    snprintf(key2,3,"d%d",dim2%10u);
    if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);
    if (sf_histfloat(in,key2,&f)) sf_putfloat(out,key1,f);

    snprintf(key1,3,"o%d",dim1);
    snprintf(key2,3,"o%d",dim2%10u);
    if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);
    if (sf_histfloat(in,key2,&f)) sf_putfloat(out,key1,f);

    snprintf(key1,7,"label%d",dim1);
    snprintf(key2,7,"label%d",dim2%10u);

    if (NULL != (val= sf_histstring(in,key1))) sf_putstring(out,key2,val);
    if (NULL != (val= sf_histstring(in,key2))) sf_putstring(out,key1,val);

    snprintf(key1,6,"unit%d",dim1);
    snprintf(key2,6,"unit%d",dim2%10u);

    if (NULL != (val= sf_histstring(in,key1))) sf_putstring(out,key2,val);
    if (NULL != (val= sf_histstring(in,key2))) sf_putstring(out,key1,val);

    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);
    
    n1=sf_esize(in);
    n2=n3=1;
    for (i=0; i < dim; i++) {
	if (i < dim1-1) {
	    n1 *= n[i]; /* block size */
	} else if (i >= dim2) {
	    n3 *= n[i]; /* loop over */
	} else {
	    n2 *= n[i]; /* read n2 blocks at a time */
	}
    }
    n12 = n1*n2;

    if (n12 + n2*sizeof(off_t) < memsize) { /* keep map incore */
	map = (off_t*) sf_alloc (n2, sizeof(off_t));
	nbuf = 0;

	for (i2=0; i2 < n2; i2++) {
	    map[i2] = make_map (dim1, dim2, n, i2);
	}

	mapfile = NULL;
    } else { /* put map out of core */
	nbuf = BUFSIZ/sizeof(off_t);
	map = (off_t*)sf_alloc(nbuf, sizeof(off_t));
	mapfile = sf_tempfile(&mapf,"w+b");

	for (i2=0, nsiz=n2; nsiz > 0; nsiz -= nbuf) {
	    if (nbuf > nsiz) nbuf=nsiz;
	    for (i=0; i < nbuf; i++, i2++) {
		map[i] = make_map (dim1, dim2, n, i2);
	    }
	    if (nbuf != fwrite(map,sizeof(off_t),nbuf,mapfile)) 
		sf_error("map write error:");
	}
    }

    if (n12 < memsize) {
	dat1 = sf_charalloc2 (n1,n2);
	dat2 = sf_charalloc2 (n1,n2);
	
	for (i3=0; i3 < n3; i3++) {
	    sf_charread(dat1[0],n12,in);
	    if (NULL == mapfile) {
		for (i2=0; i2 < n2; i2++) {
		    memcpy(dat2[i2],dat1[map[i2]],n1); 
		}
	    } else {
		if (0 > fseeko(mapfile,0,SEEK_SET))
		    sf_error ("map seek error:");

		for (i2=0, nsiz=n2; nsiz > 0; nsiz -= nbuf, i2 += nbuf) {
		    if (nbuf > nsiz) nbuf=nsiz;
		    if (nbuf != fread(map,sizeof(off_t),nbuf,mapfile)) 
			sf_error("map read error:");
		    for (i=0; i < nbuf; i++) {
			memcpy(dat2[i2+i],dat1[map[i]],n1); 
		    }
		}
	    }
	    sf_charwrite(dat2[0],n12,out);
	}
    } else {
	sf_warning("Going out of core... "
		   "(increase memsize=%zu for in-core)",memsize/(1 << 20));
	sf_unpipe(in,n12*(off_t)n3);

	buf = sf_charalloc (n1);
	
	pos = sf_tell(in);
	for (i3=0; i3 < n3; i3++) {
	    if (NULL == mapfile) {
		for (i2=0; i2 < n2; i2++) {
		    sf_seek(in,pos+(map[i2]+i3*n2)*n1,SEEK_SET);
		    sf_charread (buf,n1,in);
		    sf_charwrite(buf,n1,out);
		}
	    } else {
		if (0 > fseeko(mapfile,0,SEEK_SET))
		    sf_error ("map seek error:");

		for (nsiz=n2; nsiz > 0; nsiz -= nbuf) {
		    if (nbuf > nsiz) nbuf=nsiz;
		    if (nbuf != fread(map,sizeof(off_t),nbuf,mapfile)) 
			sf_error("map read error:");
		    for (i=0; i < nbuf; i++) {
			sf_seek(in,pos+(map[i]+i3*n2)*n1,SEEK_SET);
			sf_charread (buf,n1,in);
			sf_charwrite(buf,n1,out);
		    }
		}
	    }
	}

	if (NULL != mapfile) unlink(mapf);
    }

    exit (0);
}

static off_t make_map (int dim1, int dim2, const off_t* n, off_t i2)
{
    off_t i, j, ii[SF_MAX_DIM];

    /* from line (output) to cartesian */
    ii[dim2-1] = i2%n[dim2-1];
    j = i2/n[dim2-1];
    for (i = dim1; i < dim2-1; i++) {
	ii[i] = j%n[i];
	j /= n[i];
    }
    ii[dim1-1] = j;
    
    /* back to line (input) */
    j = ii[dim2-1];
    for (i = dim2-2; i >= dim1-1; i--) {
	j = j*n[i] + ii[i];
    }

    return j;
}
