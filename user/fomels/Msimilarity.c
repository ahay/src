/* Similarity measure between two datasets 

Takes: rect1=1 rect2=1 ... 

rectN defines the size of the smoothing stencil in N-th dimension.
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

#include "divn.h"

int main(int argc, char* argv[])
{
    int dim, dim1, i, n1, i1, i2, n2, niter;
    off_t n[SF_MAX_DIM]; 
    int rect[SF_MAX_DIM];
    char key[6];	
    float done, dtwo;
    float *one, *two, *rat1, *rat2;
    sf_file in, other, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    other = sf_input("other");

    if (SF_FLOAT != sf_gettype(in) ||
        SF_FLOAT != sf_gettype(other)) sf_error("Need float input");

    dim = sf_filedims(in,n);	

    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	if (rect[i] > 1) dim1 = i+1;
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
        if (i < dim1) {
            n1 *= n[i];
        } else {
            n2 *= n[i];
        }
    }

    if (!sf_getint("niter",&niter)) niter=20;
    /* maximum number of iterations */

    divn_init(dim1, n1, n, rect, niter);
	
    one = sf_floatalloc(n1);
    two = sf_floatalloc(n1);
    rat1 = sf_floatalloc(n1);
    rat2 = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(one,n1,in);
        sf_floatread(two,n1,other);

	/* normalization */

        done = 0.;
        dtwo = 0.;
	for (i1=0; i1 < n1; i1++) {
	    done += one[i1]*one[i1];
            dtwo += two[i1]*two[i1];
	}
        dtwo = sqrtf(n1/dtwo);	
	done = sqrtf(n1/done)/(dtwo+FLT_EPSILON);

        /* first division */

        for (i1=0; i1 < n1; i1++) {
	    one[i1] *= dtwo;
	    two[i1] *= dtwo;
	}
	
        divn(one,two,rat1);

        /* second division */

        for (i1=0; i1 < n1; i1++) {
	    one[i1] *= done;
	    two[i1] *= done;
	}
	
        divn(two,one,rat2);

	/* combination */

        for (i1=0; i1 < n1; i1++) {
           done = sqrtf(fabsf(rat1[i1]*rat2[i1]));
  	   if ((rat1[i1] > 0. && rat2[i1] < 0. && -rat2[i1] >= rat1[i1]) ||
               (rat1[i1] < 0. && rat2[i1] > 0. && rat2[i1] >= -rat1[i1])) 
		  done = -done;
           done += 1.;
           done *= done;
           done *= done/16.;
           rat1[i1] = done;	
        }

        sf_floatwrite(rat1,n1,out);
    }
 
    exit(0);
}
