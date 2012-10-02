/* Interpolation for Kirchhoff-type operators (without antialiasing). */
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

#include "stretch.h"

static int nt, nd, *x;
static float t0,dt, *w;
static bool *m;

void stretch_init (
		     int n1   /* trace length */, 
		     float o1 /* trace origin */, 
		     float d1 /* trace sampling */, 
		     int n2   /* number of data samples */)
/*< initialization >*/
{
    nt = n1; 
    t0 = o1; 
    dt = d1; 
    nd = n2;
    
    x = sf_intalloc(nd);
    m = sf_boolalloc(nd);
    w = sf_floatalloc(nd);
}

void stretch_define (const float *coord /* data coordinates [nd] */)
/*< Set up interpolation >*/
{
    int id, ix;
    float rx;

    for (id = 0; id < nd; id++) {
	    rx = coord[id];          
	    rx = (rx - t0)/dt; 
	    ix = rx; 
	    
        m[id] = ((ix < 0) || (ix > nt - 2));
	    w[id] = rx - ix;	    
	    x[id] = ix;
   }
}

void stretch_lop (bool adj    /* adjoint flag */,
		    bool add    /* addition flag */,
		    int n1, int n2, /* sizes */
		    float *ord  /* data [nd] */, 
		    float *modl /* model [nt] */)
/*< apply interpolation >*/
{
    int id, i1, i2;
    float w1, w2;

    if (n1 != nd || n2 != nt) sf_error("%s: wrong sizes",__FILE__);

    sf_adjnull(adj,add,nd,nt,ord,modl);

    for (id = 0; id < nd; id++) {
	    if (m[id]) continue;
	
	    i1 = x[id]; 
	    i2 = i1 + 1;
	    
	    w2 = w[id]; 
	    w1 = 1.0f - w2;

	    if (adj) {
		   ord[id] += w2 * modl[i2] + w1 * modl[i1];
	    } else {		
		   modl[i1] += w1 * ord[id];
		   modl[i2] += w2 * ord[id];
	    } 
	}
}

void stretch_close (void)
/*< free allocated storage >*/
{
    free (x);
    free (m);
    free (w);
}
