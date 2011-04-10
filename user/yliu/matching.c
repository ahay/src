/* Linear matching operator */
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
#include "matching.h"

static float* filt;
static int nf, n1, n2;

void matching_init (float* bb, int nff1, int nff2, int nff3) 
/*< initialize with a pointer to a matrix >*/
{
    filt   = bb;
    n1     = nff1;
    n2     = nff2;
    nf     = nff3;
}

void pmatching_lop (bool adj, bool add, 
		    int nx, int ny, float* mm, float* dd) 
/*< linear operator for Python filter structure >*/
{
    int i, j, k;
    int index;
    sf_adjnull (adj,add,nx,ny,mm,dd);

    for (k=0; k < n1; k++) {
	for (j=0; j < n2; j++) {
	    index = 0;
	    if (adj) {
		mm[j*n1+k+index] += filt[(index)*n1*n2+j*n1+k]*dd[j*n1+k];
	    } else {
		dd[j*n1+k] += filt[(index)*n1*n2+j*n1+k]*mm[j*n1+k+index];
	    }
	
	    index ++;
	    for (i=1; i < (nf+1)/2; i++) {
		
		/* zero value boundary conditions */
		if (k+i < 0 || k+i >= n1 || k-i <0 || k-i >= n1) {
		    continue; 
		}
		if (adj) {
		    mm[j*n1+k+i] += filt[(index)*n1*n2+j*n1+k]*dd[j*n1+k];
		    index ++;
		    mm[j*n1+k-i] += filt[(index)*n1*n2+j*n1+k]*dd[j*n1+k];
		    index ++;

		} else {
		    dd[j*n1+k] += filt[(index)*n1*n2+j*n1+k]*mm[j*n1+k+i];
		    index ++;
		    dd[j*n1+k] += filt[(index)*n1*n2+j*n1+k]*mm[j*n1+k-i];
		    index ++;
		}
	    }
	}
    }
}

void matching_lop (bool adj, bool add,
		   int nx, int ny, float* mm, float* dd)
/*< linear operator >*/
{
    int i, j, k;

    sf_adjnull (adj,add,nx,ny,mm,dd);

    for (k=0; k < n1; k++) {
	for (j=0; j < n2; j++) {
	    for (i=-nf/2; i < (nf+1)/2; i++) {
		/* zero value boundary conditions */
		if (k+i < 0 || k+i >= n1) {
		    continue;
		}
		if (adj) {
		    mm[j*n1+k+i] += filt[(i+nf/2)*n1*n2+j*n1+k]*dd[j*n1+k];
		} else {
		    dd[j*n1+k] += filt[(i+nf/2)*n1*n2+j*n1+k]*mm[j*n1+k+i];
		}
	    }
	}
    }
}

/* 	$Id$	 */
