/* Printing out a helix filter */
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

#include <stdlib.h>
#include <stdio.h>

#include <rsf.h>

#include "printfilter.h"
#include "boxfilter.h"

void print (int dim            /* number of dimensions */, 
	    const int *nd      /* data size [dim] */, 
	    const int *center  /* filter center [dim] */, 
	    const int *na      /* filter size [dim] */, 
	    const sf_filter aa /* filter to print */) 
/*< print a filter >*/ 
{
    float* filt; 
    int ii[SF_MAX_DIM], ii0[SF_MAX_DIM], ia, i, j, k, nf, br; 

    nf = 1;
    for (j=0; j < dim; j++) {
	nf *= na[j];
    }
    filt = sf_floatalloc(nf);

    box (dim, nd, center, na, aa, nf, filt); 
    fprintf(stderr,"-----------------------------------------\n");

    for (j=0; j < dim; j++) {
	ii0[j] = 0;
    }

    for (ia = 0; ia < nf - na[0]+1; ia += na[0]) {
	sf_line2cart(dim, na, ia, ii);
	br = 0;
	for (j=0; j < dim; j++) {
	    if (ii[j] != ii0[j]) br++;
	    ii0[j] = ii[j];
	}

        for (i=0; i < br-1; i++)
	    fprintf (stderr,"+++++++++++++++++++++++++++++++++++++++++\n");

	for (k=0; k < na[0]; k++) {
	    fprintf(stderr,"%7.3f",filt[k+ia]);
	}
	fprintf (stderr,"\n");
    }
    
    fprintf (stderr,"-----------------------------------------\n");

    free (filt);
}

/* 	$Id: printfilter.c 7107 2011-04-10 02:04:14Z ivlad $	 */
