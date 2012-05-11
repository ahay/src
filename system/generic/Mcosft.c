/* Multi-dimensional cosine transform.

The input and output are real and have the same dimensions. 
Pad the data if you need to suppress wrap-around effects.
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

#include <rsf.h>

int main (int argc, char* argv[]) 
{
    int dim, dim1, i, j, sign[SF_MAX_DIM], s[SF_MAX_DIM], n[SF_MAX_DIM];
    int n1, n2, i2, i0;
    char key[20], key2[20], *label;
    float *data, o[SF_MAX_DIM], d[SF_MAX_DIM];
    sf_file in=NULL, out=NULL;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"sign%d",i+1);
	if (1==n[i] || !sf_getint(key,sign+i)) sign[i]=0;
	/*( sign#=0 transform along #-th dimension 
	  [+1 forward or -1 backward] )*/ 
	if (sign[i]) {
	    dim1 = i;
	    snprintf(key,3,"o%d",i+1);
	    snprintf(key2,3,"u%d",i+1);
	    if (sign[i] > 0) {
		if (!sf_histfloat(in,key,o+i)) o[i]=0.;
		sf_putfloat(out,key,0.);
		sf_putfloat(out,key2,o[i]);
	    } else {		
		if (!sf_histfloat(in,key2,o+i) && 
		    !sf_getfloat(key,o+i)) o[i]=0.;
		sf_putfloat(out,key,o[i]);
	    }
	    snprintf(key,3,"d%d",i+1);
	    if (!sf_histfloat(in,key,d+i)) d[i]=1.;
	    sf_putfloat(out,key,
			1./(2*kiss_fft_next_fast_size(n[i]-1)*d[i]));

	    /* fix label and unit */
	    snprintf(key,15,"label%d",i+1);
	    snprintf(key2,15,"cosft_label%d",i+1);
	    if (NULL != (label = sf_histstring(in,key2))) {
		sf_putstring(out,key,label);
	    } else if (NULL != (label = sf_histstring(in,key))) {
		sf_putstring(out,key2,label);
		(void) sf_fft_label(i+1,label,out);
	    }
	    snprintf(key,15,"unit%d",i+1);
	    sf_fft_unit(i+1,sf_histstring(in,key),out);
	}
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    s[i] = n1;
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

    data = sf_floatalloc (n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	
	for (i=0; i <= dim1; i++) {
	    if (!sign[i]) continue;
	    sf_cosft_init(n[i] /* ,o[i],d[i] */);

	    for (j=0; j < n1/n[i]; j++) {
		i0 = sf_first_index (i,j,dim1+1,n,s);
		if (sign[i] > 0) {
		    sf_cosft_frw (data, i0, s[i]);
		} else {
		    sf_cosft_inv (data, i0, s[i]);
		}
	    }
	    sf_cosft_close();
	}
	
	sf_floatwrite(data,n1,out);
    }


    exit (0);
}
