/* Burg's method for 1-D PEF estimation */
/*
  Copyright (C) 2022 The University of Texas at Austin

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

static int nd, na;
static float *em, *s;

void burg_init(int nd1 /* data size */, 
	       int na1 /* filter size */)
/*< allocate space >*/
{
    nd = nd1;
    na = na1;

    em = sf_floatalloc(nd);
    s  = sf_floatalloc(nd);
}

void burg_close(void)
/*< free space >*/
{
    free(em);
    free(s);
}

void burg(const float *x  /* input [nd] */, 
	  float *ep       /* residual [nd] */,
	  float *a        /* filter [na] */)
/*< estimate a 1-D PEF >*/
{
    int i, ia, id;
    float c, top, bot, epi;

    for (id=0; id < nd; id++) {
	ep[id] = x[id];
	em[id] = x[id];
	s[id] = 0.0f;
    }

    a[0] = 1.0f;
    for (ia=1; ia < na; ia++) {
	a[ia] = 0.0f;
    }

    for (ia=1; ia < na; ia++) {
	top=0.0f;
	bot=0.0f;
	for (id=ia; id < nd; id++) { 
	    bot += ep[id]*ep[id]+em[id-ia]*em[id-ia];
	    top += ep[id]*em[id-ia];
	}
    
	c=2*top/(bot+SF_EPS);

	for (id=ia; id < nd; id++) {
	    epi = ep[id];
	    ep[id] -= c*em[id-ia];
	    em[id-ia] -= c*epi;
	}

	a[ia]=0.0;
	for (i=0; i < ia; i++) {
	    s[i]=a[i]-c*a[ia-i];
	}
	for (i=0; i < ia; i++) {
	    a[i]=s[i];
	}
    }
}
