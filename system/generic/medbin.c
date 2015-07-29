/* 2-D median binning */
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

static int nm, nd, *ct, *pt, *k;
static bool *m;
static float *bd;

void medbin_init (float **coord       /* coordinates [nd][2] */, 
		  float o1, float o2, 
		  float d1, float d2,
		  int   n1, int   n2  /* axes */,
		  int nd_in           /* number of data points */) 
/*< initialize >*/
{
    int id, im, i1, i2, start;

    nd = nd_in; 
    nm = n1*n2;

    ct = sf_intalloc(nm);
    pt = sf_intalloc(nm);
    k = sf_intalloc(nd);
    m = sf_boolalloc(nd);
    bd = sf_floatalloc(nd);

    for (im=0; im < nm; im++) {
	ct[im] = 0.;
    }
    for (id = 0; id < nd; id++) {
	i1 = 0.5 + (coord[id][0] - o1)/d1;
	i2 = 0.5 + (coord[id][1] - o2)/d2;
	m[id] = (bool)
	    ((i1 >= 0) && (i1 < n1) &&    
	     (i2 >= 0) && (i2 < n2));
	if (m[id]) {
	    im = i2*n1+i1;
	    ct[im]++;
	    k[id] = im; 
	}
    }

    start = 0;
    for (im=0; im < nm; im++) {
	pt[im] = start;
	start += ct[im];
    }
}

void medbin_close(void)
/*< free allocated storage >*/
{
    free(ct);
    free(pt);
    free(k);
    free(m);
    free(bd);
}

void medbin (const float *ord /* data */, 
	     float *model     /* binned */)
/*< apply median binning >*/
{
    int im, id, p, np;

    for (id = 0; id < nd; id++) { 
	if (m[id]) {
	    im = k[id];
	    bd[pt[im]] = ord[id]; 
	    pt[im]++;
	}
    }

    for (im=0; im < nm; im++) {
	np = ct[im];
	if (np > 0) {
	    p = pt[im]-np;
	    model[im] = 
		(sf_quantile(np/2,    np,bd+p)+
		 sf_quantile((np-1)/2,np,bd+p))*0.5;
	} else {
	    model[im] = 0.;
	}
    }
}
