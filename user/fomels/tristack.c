/* Resampling with triangle weights */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include "tristack.h"

static float *tmp, *tmp2, wt;
static int nc, nb, nd, np;

void tristack_init (int nbox /* triangle length */, 
		    int ndat /* coarse data length */)
/*< initialize >*/
{
    nc = ndat;
    nb = nbox;

    nd = (nc-1)*nb+1;
    np = nd + 2*nb;
    
    tmp = sf_floatalloc(np);
    tmp2 = sf_floatalloc(nc+2);

    wt = 1.0/(nb*nb);
}

void  tristack_close(void)
/*< free allocated storage >*/
{
    free (tmp);
    free (tmp2);
}

static void fold (const float *dense)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nd; i++) 
	tmp[i+nb] = dense[i];
    
    /* reflections from the right side */
    for (j=nb+nd; j < np; j += nd) {
	for (i=0; i < nd && i < np-j; i++)
	    tmp[j+i] = dense[nd-1-i];
	j += nd;
	for (i=0; i < nd && i < np-j; i++)
	    tmp[j+i] = dense[i];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nd) {
	for (i=0; i < nd && i < j; i++)
	    tmp[j-1-i] = dense[i];
	j -= nd;
	for (i=0; i < nd && i < j; i++)
	    tmp[j-1-i] = dense[nd-1-i];
    }
}

static void fold2 (float *dense)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nd; i++) 
	dense[i] += tmp[i+nb];

    /* reflections from the right side */
    for (j=nb+nd; j < np; j += nd) {
	for (i=0; i < nd && i < np-j; i++)
	    dense[nd-1-i] += tmp[j+i];
	j += nd;
	for (i=0; i < nd && i < np-j; i++)
	    dense[i] += tmp[j+i];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nd) {
	for (i=0; i < nd && i < j; i++)
	    dense[i] += tmp[j-1-i];
	j -= nd;
	for (i=0; i < nd && i < j; i++)
	    dense[nd-1-i] += tmp[j-1-i];
    }
}
    
static void doubint (void)
{
    int i, j;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=np-1; i >= 0; i--) {
	t += tmp[i];
	tmp[i] = t;
    }

    /* integrate forward */
    t=0.;
    j=0;
    for (i=0; i < np; i++) {
	t += tmp[i];
	if (0 == i%nb) {
	    tmp2[j] = t;
	    j++;
	}
    }
}

static void doubint2 (void)
{
    int i, j;
    float t;

    /* integrate backward */
    t = 0.;
    j=nc+1;
    for (i=np-1; i >= 0; i--) {
	if (0 == i%nb) {
	    t += tmp2[j];
	    j--;
	}
	tmp[i] = t;
    }

    /* integrate forward */
    t=0.;
    for (j=i=0; i < np; i++) {
	t += tmp[i];
	tmp[i] = t;
    }
}

static void triple (float* coarse)
{
    int i;
    
    for (i=0; i < nc; i++) {
	coarse[i] += (2*tmp2[i+1] - tmp2[i] - tmp2[i+2])*wt;
    }
}

static void triple2 (const float* coarse)
{
    int i;
    
    for (i=0; i < nc + 2; i++) {
	tmp2[i] = 0.;
    }

    for (i=0; i < nc; i++) {
	tmp2[i] -= coarse[i]*wt;
	tmp2[i+1] += 2*coarse[i]*wt;
	tmp2[i+2] -= coarse[i]*wt;
    }
}

void tristack (bool adj, bool add, int nx, int ny, float *dense, float *coarse)  
/*< linear operator >*/
{
    if (nd != nx || nc != ny) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj,add,nd,nc,dense,coarse);

    if (adj) {
	triple2 (coarse);
	doubint2 ();
	fold2 (dense);
    } else {
	fold (dense);
	doubint ();
	triple (coarse);
    }
}



/* 	$Id$	 */
