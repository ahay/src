/* 3D wave modeling */

/*
  Copyright (C) 2012 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#include "fd3.h"


static int nx, ny, nz, nxyz;
static int is;
static float *vel, *owv, *ud;

void wavmod_init(int n1, int n2, int n3, float dt, sf_file hvel,
	int zs, int xs, int ys)
/*< initialize >*/
{
	int i1;
	float dt2;

	nz = n1;
	nx = n2;
	ny = n3;
	nxyz = n1*n2*n3;
	is = (ys*nx+xs)*nz+zs;

	vel = sf_floatalloc(nxyz);
	owv = sf_floatalloc(nxyz);
	ud  = sf_floatalloc(nxyz);

	sf_floatread(vel, nxyz, hvel);

	dt2 = dt*dt;
	for(i1=0; i1<nxyz; i1++)
	{
		vel[i1] = vel[i1]*vel[i1]*dt2;
		owv[i1] = 0.0;
		ud[i1] = 0.0;
	}
}



void wavmod(float *wav, float source)
/*< time continuation >*/
{
	int i1;
	float t1;

	fd3_laplacian(nz, nx, ny, wav, ud);
	ud[is] += source;

	for(i1=0; i1<nxyz; i1++)
	{
		ud[i1] *= vel[i1];
		t1 = 2.0*wav[i1] - owv[i1] + ud[i1];
		owv[i1] = wav[i1];
		wav[i1] = t1;
	}
}

void wavmod_close()
/*< release memory >*/
{
	free(owv);
	free(vel);
	free(ud);
}


