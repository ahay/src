/* Create non-stationary helix filter */
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

#include "nhelix.h"
/*^*/

#include "createnhelix.h"
#include "createhelix.h"
#include "nbound.h"

nfilter createnhelix(int dim     /* number of dimensions */, 
		     int *nd     /* data size [dim] */, 
		     int *center /* filter center [dim] */, 
		     int *gap    /* filter gap [dim] */, 
		     int *na     /* filter size [dim] */, 
		     int *pch    /* patching [product(nd)] */)
/*< allocate and output a non-stationary filter >*/ 
{
    nfilter nsaa;
    sf_filter aa;
    int n123, np, ip, *nh, i;

    aa = createhelix(dim, nd, center, gap, na);

    n123=1;
    for (i=0; i < dim; i++) {
	n123 *= nd[i];
    }
    np = pch[0];
    for (i=0; i < n123; i++) {
	if (pch[i] > np) np=pch[i];
    }
    np++;

    nh = sf_intalloc(np);
    for (ip=0; ip < np; ip++) {
	nh[ip] = aa->nh;
    }
    nsaa = nallocate(np, n123, nh, pch);
    for (ip=0; ip < np; ip++) {
	for (i=0; i < aa->nh; i++) {
	    nsaa->hlx[ip]->lag[i] = aa->lag[i];
	}
	nbound(ip, dim, nd, na, nsaa); 
    }

    sf_deallocatehelix(aa);

    return nsaa;
}

