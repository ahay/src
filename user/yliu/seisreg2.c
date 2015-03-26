/* Seislet transform plus interpolation */
/*
  Copyright (C) 2008 University of Texas at Austin
   
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
#include <rsfpwd.h>

#include "seisreg2.h"
 
static int nm, nt, nd, ntm;
static float *mod, *dat, *temp;

void seisreg_init (float* coord               /* cooordinates [nd] */, 
		   float o1, float d1, int n1 /* axis */,
		   int n2                     /* number of slices */,
		   sf_interpolator interp     /* interpolation function */, 
		   int nf                     /* interpolator length */, 
		   int nd_in                  /* number of data points */,
		   bool inv                   /* inversion flag */, 
		   bool unit                  /* weighting flag */,
		   float eps                  /* regularization parameter */,
		   char type                  /* transform type */,
		   float** dd                 /* initial dip */)
/*< initialize >*/
{
    nm = n1;
    nt = n2;
    nd = nd_in;
    ntm = n1*n2;

    sf_int1_init(coord,o1,d1,n1,interp,nf,nd,0.0);
    seislet_init(n2,n1,inv,unit,eps,1,type);
    seislet_set(dd);

    mod = sf_floatalloc(nm);
    dat = sf_floatalloc(nd);
    temp = sf_floatalloc(ntm);
}

void seisreg_lop(bool adj, bool add, int nx, int ny, 
		 float *x, float *y)
/*< linear operator >*/
{
    int id, im, it;
    
    if (ny != nd*nt) sf_error("%s: wrong data size: %d != %d",__FILE__,ny,nd*nt);
    if (nx != nm*nt) sf_error("%s: wrong data size: %d != %d",__FILE__,nx,nm*nt);

    sf_adjnull (adj,add,nx,ny,x,y);

    if (adj) {
	for (it=0; it < nt; it++) { /* loop over time slices */
	    for (id=0; id < nd; id++) {
		dat[id] = y[id*nt+it];
	    }
	    /* apply interpolation */
	    sf_int1_lop(adj,false,nm,nd,mod,dat);
	    for (im=0; im < nm; im++) {
		temp[im*nt+it] = mod[im];
	    }
	}
	seislet_construct(adj,false,ntm,ntm,x,temp);

    } else {
	seislet_construct(adj,false,ntm,ntm,x,temp);
	for(it=0; it < nt; it++) { /* loop over time slices */
	    for (im=0; im < nm; im++) {
		mod[im] = temp[im*nt+it];
	    }
	    /* apply interpolation */
	    sf_int1_lop(adj,false,nm,nd,mod,dat);
	    for (id=0; id < nd; id++) {
		y[id*nt+it] = dat[id];
	    }
	} 
    }
}

void seisreg_close(void)
/*< deallocate space >*/
{
    free(mod);
    free(dat);
    free(temp);
    sf_int1_close();
    seislet_close();
}

/* 	$Id$	 */
