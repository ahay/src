/* Convert input to its derivative using LS and shaping regularization
 * applied to causint_lop */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
    int i, niter, nd, dim, id;
    int n[SF_MAX_DIM], box[SF_MAX_DIM];
    float *d, *m, *wt, *m0, wti;
    char key[6];
    sf_file data, modl, weight;

    sf_init(argc,argv);
    data = sf_input("in");
    modl = sf_output("out");
    
    if (NULL != sf_getstring("weight")) {
	weight = sf_input("weight");
    } else {
	weight = NULL;
    }

    dim = sf_filedims (data,n);

    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
    }

    for (i=0; i < dim; i++) { 	 
	snprintf(key,6,"rect%d",i+1); 	 	 
	if (!sf_getint(key,box+i)) box[i]=1; 	 
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
    } 	 
	 
    smoothder_init(nd, dim, box, n);

    d = sf_floatalloc(nd);
    m = sf_floatalloc(nd);
    wt = sf_floatalloc(nd);
    m0 = sf_floatalloc(nd);

    sf_floatread(d,nd,data);

    if (NULL != weight) {
	sf_floatread(wt,nd,weight);

	wti = 0.;
	for (id=0; id < nd; id++) {
	    wti += wt[id]*wt[id];
	}
	if (wti > 0.) wti = sqrtf(nd/wti);
	for (id=0; id < nd; id++) {
	    wt[id] *= wti;
	}
    } else {
	for (id=0; id < nd; id++) {
	    wt[id] = 1.0;
	}	
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    for (id=0; id < nd; id++) {
	m0[id] = -d[id];
    }
    
    sf_repeat_lop(false,true,nd,nd,m0,d);

    smoothder(niter, wt, d, m);
 
    for (id=0; id < nd; id++) {
	m[id] -= m0[id];
    }

    sf_floatwrite(m,nd,modl);

    exit(0);
}

