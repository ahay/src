/* Linearized eikonal operator */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

#include "lineiko.h"

static int n, *sq, ***cd;
static float d[3], *sref, *tref, *timex=NULL, *slow;

void lineiko_init (int n12      /* total number of grid points */, 
		   int *sq1     /* sequence array */, 
		   int ***cd1   /* dependence matrix */, 
		   float *d1    /* [3] grid sampling */, 
		   float *sref1 /* reference slowness */, 
		   float *tref1 /* reference traveltime */)
/*< initialize >*/
{
    n = n12; 
    d[0] = 1./(d1[0]*d1[0]);
    d[1] = 1./(d1[1]*d1[1]);
    d[2] = 1./(d1[2]*d1[2]);

    sq = sq1; 
    cd = cd1; 
    sref = sref1; 
    tref = tref1;

    if (NULL == timex) {
	timex = sf_floatalloc(n12);
	slow  = sf_floatalloc(n12);
    }
}

void lineiko_close(void)
/*< free allocated storage >*/
{
    if (NULL != timex) {
	free(timex);
	free(slow);
	timex = NULL;
    }
}

void lineiko_lop (bool adj, bool add, int ns, int nt, 
		  float *ss /* slowness perturbation */, 
		  float *tt /* traveltime pertuurbation */)
/*< linear operator >*/
{
    int iq, i, j, k, *from;
    float a, b, s, s2, sten[3];

    if (n != ns || n != nt) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull (adj,add,ns,nt,ss,tt);
    
    if (adj) {
	for (j=0; j < n; j++) {
	    slow[j] = 0.;
	}

	for (iq = n-1; iq >=0; iq--) {
	    i = sq[iq];
	    from = cd[i][0];

	    a = 0.;
	    for (k=0; k < 3; k++) {
		j = from[k];
		if (j >= 0) {
		    sten[k] = tref[i] - tref[j];
		    a += sten[k]*d[k];
		}
	    }

	    slow[i] = sref[i]*(slow[i] + tt[i])*a/(a*a + FLT_EPSILON);

	    for (k=0; k < 3; k++) {
		j = from[k];
		if (j >= 0) {
		    slow[j] += slow[i]*sten[k]*d[k]/sref[i];
		}
	    }
	}

	for (j=0; j < n; j++) {
	    ss[j] += slow[j];
	}
    } else {
	for (j=0; j < n; j++) {
	    timex[j] = 0.;
	}

	for (iq=0; iq < n; iq++) {
	    i = sq[iq];
	    from = cd[i][0];

	    a = 0.;
	    b = ss[i]*sref[i];
	    for (k=0; k < 3; k++) {
		j = from[k];
		if (j >= 0) {
		    s = timex[j];
		    s2 = tref[i] - tref[j];
		    
		    a += s2*d[k];
		    b += s*s2*d[k];
		}
	    }
	              
	    timex[i] = b*a/(a*a + FLT_EPSILON);
	}

	for (j=0; j < n; j++) {
	    tt[j] += timex[j];
	}
    }
}

