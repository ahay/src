/* Simple reflection tomography  */
/*
 Copyright (C) 2008 Colorado School of Mines
 
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

#include "xtomo.h"

static int nt,nz,nx,nh,ny;
static float oz,dz,ox,dx,oh,dh,oy,dy;

void xtomo_init(float oz1, float dz1, 
		float ox1, float dx1, 
		float oh1, float dh1,
		float oy1, float dy1,
		int nt1, int nz1, int nx1, int nh1, int ny1)
/*< initialize with domensions >*/
{
    oz = oz1; ox = ox1; oh = oh1; oy = oy1; 
    dz = dz1; dx = dx1; dh = dh1; dy = dy1;
    nz = nz1; nx = nx1; nh = nh1; ny = ny1;
} 

void xtomo_lop(bool adj, bool add, int nm, int nd, float *modl, float *data)
/*< linear operator >*/
{
    float z, x, h, y, zmax;
    int iz, ix, ih, iy, it, id, im;

    if (nm != nt*nz*nx || nd != nt*nh*ny) 
	sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj, add, nm, nd, modl, data);

    zmax=oz+dz*(nz-1);
    for (ix=0; ix < nx; ix++) {
	x = ox + ix*dx;
	for (ih=0; ih < nh; ih++) {
	    h = oh + ih*dh;
	    for (iz=0; iz < nz; iz++) {
		z = oz + iz*dz;

		y = x +  h*z/zmax; /* maybe opposite sign on h or (1-z/zmax) */

		iy = 0.5 + (y-oy)/dy;

		if(0 <= iy && iy <= ny-1) {
		    for (it=0; it < nt; it++) {
			id = it+nt*(ih+iy*nh);
			im = it+nt*(iz+ix*nz);

			if (adj) {
			    modl[im] += data[id];
			} else {
			    data[id] += modl[im];
			}
		    }
		}
	    }
	}
    }
}
