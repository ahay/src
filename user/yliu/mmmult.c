/* Multi-dimension matrix multiplication operator */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
/*^*/
#include "mmmult.h"

static float**** filt;
static int nf1, nf2, n1, n2;

void matmult_init (float**** bb, int nff1, int nff2, int nff3, int nff4) 
/*< initialize with a pointer to a matrix >*/
{
    filt = bb;
    nf1  = nff1;
    nf2  = nff2;
    n1   = nff3;
    n2   = nff4;
}

void matmult_lop (bool adj, bool add, 
		  int nx, int ny, float** model, float** data) 
/*< linear operator >*/
{
    int i, j, k, l;
    
    sf_adjnull (adj,add,nx,ny,model[0],data[0]);
    
    for (l=0; l < n2; l++) {
	for (k=0; k < n1; k++) {
	    for (j=0; j < nf2; j++) {
		for (i=-nf1/2; i < (nf1-nf1/2-1); i++) {
		    /* zero value boundary conditions */
		    if (l+j < 0 || l+j >= n2 || k+i < 0 || k+i >= n1) {
			continue; 
		    }
		    if (adj) {
			model[l+j][k+i] += filt[l][k][j][i]*data[l][k];
		    } else {
			data[l][k] += filt[l][k][j][i]*model[l+j][k+i];
		    }
		}
	    }
	}
    }
}

/* 	$Id$	 */
