/* Window a portion of a dataset. 

Other parameters from the command line are passed to the output (similar to sfput).
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

static off_t seektable (int dim, off_t *n, off_t *m, off_t *f, off_t *j, 
			off_t n1, off_t n2, off_t *table);

int main (int argc, char *argv[])
{
    int i, esize, dim, tmp;
    off_t n1, n2, m1, i2, i1, j1, jump;
    off_t i0, j[SF_MAX_DIM];
    off_t n[SF_MAX_DIM], m[SF_MAX_DIM], f[SF_MAX_DIM], *table, maxsize;
    float a, d[SF_MAX_DIM], o[SF_MAX_DIM];
    char key[8], *label[SF_MAX_DIM], *unit[SF_MAX_DIM], *buf;
    bool squeeze, verb;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint (in,"esize",&esize) || esize <= 0)
	sf_error("Need esize > 0 in in");

    dim = sf_largefiledims(in,n);

    for (i=0; i < dim; i++) {
	/* get o's */
	snprintf(key,4,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	
        /* get d's */
	snprintf(key,4,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) d[i]=1.;

	/* get j's */
	snprintf(key,4,"j%d",i+1);
	if (!sf_getint(key,&tmp)) {
	    /*( j#=(1,...) jump in #-th dimension )*/
	    snprintf(key,4,"d%d",i+1);
	    if (sf_getfloat(key,&a)) {
		/*( d#=(d1,d2,...) sampling in #-th dimension )*/
		j[i] = 0.5 + a/d[i];
	    } else {
		j[i] = 1;
	    }
	} else
            j[i] = tmp;

	/* get f's */	
	snprintf(key,4,"f%d",i+1);
	if (!sf_getlargeint(key,f+i)) {
	    /*( f#=(0,...) window start in #-th dimension )*/
	    snprintf(key,6,"min%d",i+1);
	    if (sf_getfloat(key,&a)) {
		/*( min#=(o1,o2,,...) minimum in #-th dimension )*/
		f[i] = 0.5 + (a - o[i]) / d[i];
	    } else {
		f[i] = 0;
	    }
	}
	if (f[i] < 0) {
	    f[i] = n[i]+f[i];
	    if (f[i] < 0) sf_error("Negative f%d=%d",i+1,f[i]);
	}

	/* new values for o and d */
	o[i] += f[i]*d[i]; 	
	d[i] *= j[i];

	/* get n's */
	snprintf(key,4,"n%d",i+1);
	if (!sf_getlargeint(key,m+i)) {
	    /*( n#=(0,...) window size in #-th dimension )*/
	    snprintf(key,6,"max%d",i+1);
	    if (sf_getfloat(key,&a)) {
		/*( max#=(o1+(n1-1)*d1,o2+(n1-1)*d2,,...) maximum in #-th dimension )*/
		m[i] = 1.5 + (a - o[i]) / d[i];
	    } else {
		m[i] = 1.5 + (n[i] - 1 - f[i]) / j[i];
	    }
	}
	if (f[i]+(m[i]-1)*j[i] > n[i]) 
#if defined(__cplusplus) || defined(c_plusplus)
	    sf_error ("m%d=%ld is too big",i+1,(long) m[i]);
#else
            sf_error ("m%d=%lld is too big",i+1,(long long) m[i]);	
#endif

	/* get labels */
	snprintf(key,8,"label%d",i+1);
	label[i] = sf_histstring(in,key);

	/* get units */
	snprintf(key,7,"unit%d",i+1);
	unit[i] = sf_histstring(in,key);
    }

    if (!sf_getbool("verb",&verb)) verb=false;
    /* Verbosity flag */

    if (verb) {
	for (i=0; i < dim; i++) {
	    if (m[i] != n[i])
 #if defined(__cplusplus) || defined(c_plusplus)
		sf_warning("Windowing f%d=%ld j%d=%d n%d=%ld"
			   " min%d=%g max%d=%g",
			   i+1,(long) f[i],i+1,j[i],i+1,(long) m[i],
			   i+1,o[i],i+1,o[i]+(m[i]-1)*d[i]);
#else
	        sf_warning("Windowing f%d=%lld j%d=%d n%d=%lld"
			   " min%d=%g max%d=%g",
			   i+1,(long long) f[i],i+1,j[i],i+1,(long long) m[i],
			   i+1,o[i],i+1,o[i]+(m[i]-1)*d[i]);
#endif
	}
    }

    if (!sf_getbool("squeeze",&squeeze)) squeeze=true;
    /* if y, squeeze dimensions equal to 1 to the end */

    sf_expandpars(out);

    for (i=i0=0; i0 < dim; i0++) {
	if (squeeze && 1==m[i0]) continue;
	snprintf(key,4,"n%d",(i+1)%10u);
	sf_putlargeint(out,key,m[i0]);
	snprintf(key,4,"o%d",(i+1)%10u);
	sf_putfloat(out,key,o[i0]);
	snprintf(key,4,"d%d",(i+1)%10u);
	sf_putfloat(out,key,d[i0]);
	if (NULL != label[i0]) {
	    snprintf(key,8,"label%d",(i+1)%10u);
	    sf_putstring(out,key,label[i0]);
	}
	if (NULL != unit[i0]) {
	    snprintf(key,7,"unit%d",(i+1)%10u);
	    sf_putstring(out,key,unit[i0]);
	}
	i++;
    }

    if (squeeze) {
	for (i0=0; i0 < dim; i0++) {
	    if (1 != m[i0]) continue;
	    snprintf(key,4,"n%d",(i+1)%10u);
	    sf_putlargeint(out,key,m[i0]);
	    snprintf(key,4,"o%d",(i+1)%10u);
	    sf_putfloat(out,key,o[i0]);
	    snprintf(key,4,"d%d",(i+1)%10u);
	    sf_putfloat(out,key,d[i0]);
	    if (NULL != label[i0]) {
		snprintf(key,8,"label%d",(i+1)%10u);
		sf_putstring(out,key,label[i0]);
	    }
	    if (NULL != unit[i0]) {
		snprintf(key,7,"unit%d",(i+1)%10u);
		sf_putstring(out,key,unit[i0]);
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

    buf = sf_charalloc (n1);
    table = sf_largeintalloc (n2);

    maxsize = seektable(dim,n,m,f,j,n1,n2,table);

    if (verb) sf_warning("maxsize=%zu",maxsize);

    sf_unpipe(in,maxsize);

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


    exit (0);
}

static off_t seektable(int dim, off_t *n, off_t *m, off_t *f, off_t *j, 
		       off_t n1, off_t n2, off_t *table)
{
    off_t i2, i, ii[SF_MAX_DIM];
    off_t t, t2;

    t2 = sf_large_cart2line (dim-1, n+1, f+1);
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

    return (t2*n[0]+f[0]+n1);
}
