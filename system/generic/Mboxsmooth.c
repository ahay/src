/* Multi-dimensional smoothing with boxes. */
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
    int dim, dim1, i, j, k, rect[SF_MAX_DIM], s[SF_MAX_DIM], n[SF_MAX_DIM];
    int nrep, irep, n1, n2, i2, i0, rmax;
    char key[6];
    float *data, *data2;
    sf_file in=NULL, out=NULL;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    rmax = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
        /*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) {
	    dim1 = i;
	    if (rect[i]+n[i] > rmax) rmax = rect[i]+n[i];
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
    data2 = sf_floatalloc (rmax);

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	for (i=0; i <= dim1; i++) {
	    if (rect[i] <= 1) continue;
	    sf_box_init (rect[i],n[i],false);
	    for (j=0; j < n1/n[i]; j++) {
		i0 = sf_first_index (i,j,dim1+1,n,s);
		for (irep=0; irep < nrep; irep++) {
		    sf_boxsmooth2 (i0,s[i],data,data2);
		    for (k=0; k < n[i]; k++) {
			data[i0+k*s[i]] = data2[k];
		    }
		}
	    }
	}
	
	sf_floatwrite(data,n1,out);
    }


    exit (0);
}
