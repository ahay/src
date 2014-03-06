/* Spraying linear operator*/
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

static int ns, ns2;

int spray_init(int nr /* spray radius */)
/*< initialize >*/
{
    ns=nr;
    ns2=2*ns+1;
    return ns2;
}

void spray_lop(bool adj, bool add, int n, int nu, float* trace, float *u)
/*< linear operator >*/
{
    int i, is, ip, i1;

    if (nu != n*ns2) sf_error("%s: wrong size %d != %d*%d",__FILE__,nu,n,ns2);

    sf_adjnull(adj,add,n,nu,trace,u);

    for (i=0; i < n; i++) { 
	if (adj) {
	    trace[i] += u[i*ns2+ns];
	} else {	    
	    u[i*ns2+ns] += trace[i];
	}

	/* predict forward */
	for (is=0; is < ns; is++) {
	    ip = i-is-1;
	    if (ip < 0) break;
	    i1 = ip*ns2+ns-is-1;
	    if (adj) {
		trace[i] += u[i1];
	    } else {
		u[i1] += trace[i];
	    }
	}
	    
	/* predict backward */
	for (is=0; is < ns; is++) {
	    ip = i+is+1;
	    if (ip >= n) break;
	    i1 = ip*ns2+ns+is+1;
	    if (adj) {
		trace[i] += u[i1];
	    } else {
		u[i1] += trace[i];
	    }
	}
    }
}

