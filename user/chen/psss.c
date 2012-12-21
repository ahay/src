/* step split phase shift */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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

static int nw,nx;
static float dz;
static float *vel;
static int *k2;

void sf_psss_init(int nw0,int nx0,int nz,
	float dw,float dx, float dz0, float *v0)
/*< initialize >*/
{
	int i;
	nw=nw0;
	nx=nx0;

	dz=dz0*2.0*SF_PI/(dx*nx);

	vel=(float*)sf_floatalloc(nz);
	k2=(int*)sf_intalloc(nx);

	for(i=0; i<nx/2+1; i++)		k2[i] = i*i;
	for(i=nx/2+1; i<nx; i++) 	k2[i] = (i-nx)*(i-nx);

	for(i=0; i<nz; i++)
	{
		vel[i] = dw*dx*nx/v0[i];
		vel[i] *= vel[i];
	}
}

void sf_psss_step(int iz, sf_complex **io) 
/*< step in depth >*/
{
	int iw,ix;
	float t1;
	sf_complex a;

	for(ix=0;ix<nx;ix++)
	{
		for(iw=0;iw<nw;iw++)	
		{
			t1  = (vel[iz]*iw*iw-k2[ix]);
			if(t1>0) 
			{
				t1 = sqrt(t1)*dz;
				a = sf_cmplx(cos(t1),sin(t1));
			}else	a = 0.0;
			io[ix][iw] *= a;
		}
	}
}

void sf_psss_close()
/*< free allocated storage >*/
{
	free(vel);
}

