/* Convert input to its derivative using LS and shaping regularization
 * applied to causint_lop d = L m

 Takes: rect1= rect2= ...

 rectN defines the size of the smoothing stencil in N-th dimension.
*/
/*
  Copyright (C) 2010 Politecnico di Milano
  
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

#include <math.h>

#include <rsf.h>

#include "smoothder.h"

int main(int argc, char* argv[])
{
    int i, niter, nd, dim, n1, n2, i1, i2;
    int n[SF_MAX_DIM], box[SF_MAX_DIM];
    float **data, **model, *model0; 
    char key[6];
    sf_file DATA, MODEL, DATA_OUT;

    sf_init(argc,argv);
    DATA = sf_input("in");
    MODEL = sf_output("out");

    dim = sf_filedims (DATA,n);

    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
    }
    n1 = n[0];
    n2 = nd/n1;

    for (i=0; i < dim; i++) { 	 
	snprintf(key,6,"rect%d",i+1); 	 	 
	if (!sf_getint(key,box+i)) box[i]=1; 	 
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
    } 	 
	 
    smoothder_init(nd, dim, box, n);

    data = sf_floatalloc2(n1,n2);
    model = sf_floatalloc2(n1,n2);
    /* model0 = sf_floatalloc2(n1,n2); */
    model0 = sf_floatalloc(n2);
    /* wt = sf_floatalloc2(n1,n2); */

    sf_floatread(data[0],nd,DATA);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */


    for (i2=0; i2 < n2; i2++) {
    	model0[i2] = data[i2][0];
	for (i1=0; i1 < n1; i1++) {
	    /* model0[i2][i1] = -data[i2][0]; */
	    data[i2][i1] = data[i2][i1]-model0[i2];
/*	   fprintf(stderr,"data = %f\n n1=%d",data[i2][i1],n1); */

	}
    }
/*    fprintf(stderr,"data = %f",data[0][0]); */
    /* sf_repeat_lop(false,true,nd,nd,model0[0],data[0]); */
    smoothder(niter, NULL, data[0], model[0]);
 
    sf_floatwrite(model[0],nd,MODEL);

    if (NULL != sf_getstring("dataout")) {
	/* optionally, output predicted data */
	sf_repeat_lop(false,false,nd,nd,model[0],data[0]);

    	for (i2=0; i2<n2;i2++) {
	    for (i1=0;i1<n1;i1++) {
		data[i2][i1] = data[i2][i1]+model0[i2];
	    }
    	}
	/* predicted or recalculated data*/
    	DATA_OUT = sf_output("dataout");

	sf_floatwrite(data[0],nd,DATA_OUT);
    }

    exit(0);
}

/* 	$Id: Mdix.c 5595 2010-03-21 16:54:14Z sfomel $	 */
