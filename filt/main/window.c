/* Window a portion of the dataset. 
   
Takes: [j1=1 j2=1 ... f1=0 f2=0 ... n1=n1 n2=n2 ... max1= max2= ... min1= min2= ...]

jN defines the jump in N-th dimension
fN is the window start
nN is the window size
minN and maxN is the maximum and minimum in N-th dimension
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

#include <rsf.h>

static void seektable (int dim, int *n, int *m, int *f, int *j, 
		       int n1, int n2, off_t *table);

int main (int argc, char *argv[])
{
    int i, esize, dim, n1, n2, m1, i2, i1, j1, jump;
    int i0, n[SF_MAX_DIM], m[SF_MAX_DIM], j[SF_MAX_DIM], f[SF_MAX_DIM];
    off_t *table;
    float a, d[SF_MAX_DIM], o[SF_MAX_DIM];
    char key[7], *label[SF_MAX_DIM], *buf;
    bool squeeze, verb;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint (in,"esize",&esize) || esize <= 0)
	sf_error("Need esize > 0 in in");

    dim = sf_filedims(in,n);

    for (i=0; i < dim; i++) {
	/* get o's */
	snprintf(key,3,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	
        /* get d's */
	snprintf(key,3,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) d[i]=1.;

	/* get j's */
	snprintf(key,3,"j%d",i+1);
	if (!sf_getint(key,j+i)) {
	    snprintf(key,3,"d%d",i+1);
	    if (sf_getfloat(key,&a)) {
		j[i] = 0.5 + a/d[i];
	    } else {
		j[i] = 1;
	    }
	} 

	/* get f's */	
	snprintf(key,3,"f%d",i+1);
	if (!sf_getint(key,f+i)) {
	    snprintf(key,5,"min%d",i+1);
	    if (sf_getfloat(key,&a)) {
		f[i] = 0.5 + (a - o[i]) / d[i];
	    } else {
		f[i] = 0;
	    }
	}
	if (f[i] < 0) sf_error("Negative f%d=%d",i+1,f[i]);

	/* new values for o and d */
	o[i] += f[i]*d[i]; 	
	d[i] *= j[i];

	/* get n's */
	snprintf(key,3,"n%d",i+1);
	if (!sf_getint(key,m+i)) { 
	    snprintf(key,5,"max%d",i+1);
	    if (sf_getfloat(key,&a)) {
		m[i] = 1.5 + (a - o[i]) / d[i];
	    } else {
		m[i] = 1.5 + (n[i] - 1 - f[i]) / j[i];
	    }
	}
	if (1+(m[i]-1)*j[i] > n[i]) 
	    sf_error ("n%d=%d is too big",i+1,m[i]);

	/* get label's */
	snprintf(key,7,"label%d",i+1);
	label[i] = sf_histstring(in,key);
    }

    if (!sf_getbool("verb",&verb)) verb=false;
    /* Verbosity flag */

    if (verb) {
	for (i=0; i < dim; i++) {
	    if (m[i] != n[i]) 
		sf_warning("Windowing f%d=%d j%d=%d n%d=%d min%d=%g max%d=%g",
			   i+1,f[i],i+1,j[i],i+1,m[i],
			   i+1,o[i],i+1,o[i]+(m[i]-1)*d[i]);
	}
    }

    if (!sf_getbool("squeeze",&squeeze)) squeeze=true;
    /* if y, squeeze dimensions equal to 1 to the end */

    for (i=i0=0; i0 < dim; i0++) {
	if (squeeze && 1==m[i0]) continue;
	snprintf(key,3,"n%d",i+1);
	sf_putint(out,key,m[i0]);
	snprintf(key,3,"o%d",i+1);
	sf_putfloat(out,key,o[i0]);
	snprintf(key,3,"d%d",i+1);
	sf_putfloat(out,key,d[i0]);
	if (NULL != label[i0]) {
	    snprintf(key,7,"label%d",i+1);
	    sf_putstring(out,key,label[i0]);
	}
	i++;
    }

    if (squeeze) {
	for (i0=0; i0 < dim; i0++) {
	    if (1 != m[i0]) continue;
	    snprintf(key,3,"n%d",i+1);
	    sf_putint(out,key,m[i0]);
	    snprintf(key,3,"o%d",i+1);
	    sf_putfloat(out,key,o[i0]);
	    snprintf(key,3,"d%d",i+1);
	    sf_putfloat(out,key,d[i0]);
	    if (NULL != label[i0]) {
		snprintf(key,7,"label%d",i+1);
		sf_putstring(out,key,label[i0]);
	    }
	    i++;
	}
    }

    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);
    
    /* Now do the actual work */
    n2 = sf_filesize(out)/m[0];
    m1 = m[0]*esize;
    n1 = (1+(m[0]-1)*j[0])*esize;
    jump = (j[0]-1) * esize;
    n[0] *= esize;
    f[0] *= esize;

    sf_unpipe(in,sf_filesize(in)*esize);

    buf = sf_charalloc (n1);
    table = (off_t*) sf_alloc (n2,sizeof(off_t));

    seektable(dim,n,m,f,j,n1,n2,table);

    for (i2=0; i2 < n2; i2++) {
	if (table[i2]) sf_seek(in,table[i2],SEEK_CUR);
	sf_charread(buf,n1,in);
	if (jump) {
	    for (i1=j1=0; i1 < m1; j1 += jump) {
		for (i=0; i < esize; i++, i1++, j1++) {
		    buf[i1] = buf[j1];
		}
	    }
	}

	sf_charwrite(buf,m1,out);
    }

    sf_close();
    exit (0);
}

static void seektable(int dim, int *n, int *m, int *f, int *j, 
		      int n1, int n2, off_t *table)
{
    int i2, t, t2, i, ii[SF_MAX_DIM];

    t2 = sf_cart2line (dim-1, n+1, f+1);
    table[0] = t2*n[0] + f[0];

    for (i2=1; i2 < n2; i2++) {
	t = i2;
	for (i = 1; i < dim-1; i++) {
	    /* cartesian coordinates in window */
	    ii[i] = t%m[i];
	    t /= m[i];
	}
	t = f[dim-1] + t*j[dim-1];
	for (i = dim-2; i >= 1; i--) {
	    /* line coordinates in input */
	    t = t*n[i] + f[i] + ii[i]*j[i];
	}
	table[i2] = (t-t2)*n[0]-n1; 
	t2 = t;
    }
}

/* 	$Id$	 */
