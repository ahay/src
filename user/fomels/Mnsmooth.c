/* N-D non-stationary smoothing. 

Takes: rect1=rect1.rsf rect2=rect2.rsf ... 

rectN defines the sizes of the smoothing stencils in N-th dimension.
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

#include "ntrianglen.h"

int main (int argc, char* argv[]) 
{
    int *rct[SF_MAX_DIM], box[SF_MAX_DIM], n[SF_MAX_DIM], s[SF_MAX_DIM];
    int dim, dim1, i, n1, n2, i1, i2;
    float *data, *smoo;
    char key[6];
    sf_file in, out, rect[SF_MAX_DIM];

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
 
    dim = sf_filedims (in,n);

    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (NULL != sf_getstring(key)) {
	    rect[i] = sf_input(key);
	    if (SF_INT != sf_gettype(rect[i])) sf_error("Need int %s",key);

	    dim1 = i;
	} else {
	    rect[i] = NULL;
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
    smoo = sf_floatalloc (n1);

    for (i=0; i <= dim1; i++) {
	box[i] = 1;
	if (NULL != rect[i]) {
	    rct[i] = sf_intalloc (n1);
	    sf_intread(rct[i],n1,rect[i]);
	    sf_fileclose(rect[i]);

	    for (i1=0; i1 < n1; i1++) {
		if (rct[i][i1] > box[i]) box[i] = rct[i][i1];
	    }
	} else {
	    rct[i] = NULL;
	}
    }

    ntrianglen_init(dim1+1,box,n,rct);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	ntrianglen_lop(false,false,n1,n1,data,smoo);
	
	sf_floatwrite(smoo,n1,out);
    }    

    exit (0);
}

/* 	$Id: Msmooth.c 691 2004-07-04 19:28:08Z fomels $	 */
