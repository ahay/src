/* Transpose two axes in a dataset. */
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

static void make_map (int dim1, int dim2, 
		      const int* n, int n2, int* map);

int main(int argc, char* argv[])
{
    int i, dim, n[SF_MAX_DIM], n1, n2, n3;
    int dim1, dim2, i2, i3, *map;
    const int mem=100;
    off_t pos, memsize;
    char key1[7], key2[7], *val, **dat1, **dat2, *buf;
    sf_file in, out;
    float f;

    sf_init (argc,argv);
    in  = sf_input  ( "in");
    out = sf_output ("out");

    memsize = sf_memsize(mem);
    /* Available memory size (in Mb) */
    
    dim = sf_filedims(in,n);

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
    snprintf(key2,3,"n%d",dim2);
    sf_putint(out,key1,n[dim2-1]);
    sf_putint(out,key2,n[dim1-1]);
 
    snprintf(key1,3,"d%d",dim1);
    snprintf(key2,3,"d%d",dim2);
    if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);
    if (sf_histfloat(in,key2,&f)) sf_putfloat(out,key1,f);

    snprintf(key1,3,"o%d",dim1);
    snprintf(key2,3,"o%d",dim2);
    if (sf_histfloat(in,key1,&f)) sf_putfloat(out,key2,f);
    if (sf_histfloat(in,key2,&f)) sf_putfloat(out,key1,f);

    snprintf(key1,7,"label%d",dim1);
    snprintf(key2,7,"label%d",dim2);

    if (NULL != (val= sf_histstring(in,key1))) sf_putstring(out,key2,val);
    if (NULL != (val= sf_histstring(in,key2))) sf_putstring(out,key1,val);

    snprintf(key1,6,"unit%d",dim1);
    snprintf(key2,6,"unit%d",dim2);

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

    map = sf_intalloc (n2);
    make_map (dim1, dim2, n, n2, map);

    if ((off_t) n1*n2 < memsize) {
	dat1 = sf_charalloc2 (n1,n2);
	dat2 = sf_charalloc2 (n1,n2);
	
	for (i3=0; i3 < n3; i3++) {
	    sf_charread(dat1[0],n1*n2,in);
	    for (i2=0; i2 < n2; i2++) {
		memcpy(dat2[i2],dat1[map[i2]],n1); 
	    }
	    sf_charwrite(dat2[0],n1*n2,out);
	}
    } else {
	sf_warning("Going out of core... "
		   "(increase memsize=%d for in-core)",memsize/(1 << 20));
	sf_unpipe(in,(off_t) n1*n2*n3);

	buf = sf_charalloc (n1);
	
	pos = sf_tell(in);
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		sf_seek(in,pos+(off_t) (map[i2]+i3*n2)*n1,SEEK_SET);
		sf_charread (buf,n1,in);
		sf_charwrite(buf,n1,out);
	    }
	}
	
	sf_close();
    }

    exit (0);
}

static void make_map (int dim1, int dim2,
		      const int* n, int n2, int* map)
{
    int i, j, i2;
    int ii[SF_MAX_DIM];

    for (i2=0; i2 < n2; i2++) {

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
	map[i2] = j;
    }
}
