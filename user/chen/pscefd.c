/* explicit finite difference phase shift */

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

static int nw,nx,cas;
static float f0,dz,dx2;
static float **vel,*buf;

void sf_pscefd_init(int cascade,int nx0,int nw0,
	float dx,float dz0,
	float f00,float **v0)
/*< initialize >*/
{
	cas=cascade;
	nx=nx0;
	nw=nw0;

	dz=dz0/cas;
	f0=f00;
	vel=v0;

	dx2=(dx*dx);
	buf=sf_floatalloc(nx*2);
}


static void pscefd_step3(int iz,float **io) 
{
	int iw,ix;
	float t1,a0i,a1i,a2i, br,bi,cr,ci;

	for(iw=1;iw<nw;iw++)  /* frequency slice */
	{
		t1=dz/(iw*dx2*f0);
		/* ix=0 boundary */
		ix=0;
		a1i=-t1*vel[ix][iz];
		a2i=0.5*t1*vel[ix+1][iz];
		br=f0*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[iw][ix*2]
			- io[iw][ix*2+1] * a1i
			- io[iw][ix*2+3] * a2i;
		ci=	+ io[iw][ix*2+1]
			+ io[iw][ix*2] * a1i
			+ io[iw][ix*2+2] * a2i;
		buf[ix*2]   = cr*br - ci*bi;	
		buf[ix*2+1] = cr*bi + ci*br;		

		for(ix=1;ix<nx-1;ix++)
		{        
			a0i=0.5*t1*vel[ix-1][iz];
			a1i=-t1*vel[ix][iz];
			a2i=0.5*t1*vel[ix+1][iz];
			br=f0*dz*iw/vel[ix][iz];
			bi=sin(br);
			br=cos(br);
			cr=	  io[iw][ix*2] 
				- io[iw][ix*2-1] * a0i
				- io[iw][ix*2+1] * a1i
				- io[iw][ix*2+3] * a2i;
			ci=	+ io[iw][ix*2+1]
				+ io[iw][ix*2-2] * a0i
				+ io[iw][ix*2] * a1i
				+ io[iw][ix*2+2] * a2i;
			buf[ix*2]   = cr*br - ci*bi;	
			buf[ix*2+1] = cr*bi + ci*br;	
		}

		/* ix=nx-1 boundary */
		a0i=0.5*t1*vel[ix-1][iz];
		a1i=-t1*vel[ix][iz];
		br=f0*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[iw][ix*2]
			- io[iw][ix*2-1] * a0i
			- io[iw][ix*2+1] * a1i;
		ci=	+ io[iw][ix*2+1] 
			+ io[iw][ix*2] * a1i
			+ io[iw][ix*2+2] * a1i;
		buf[ix*2]   = cr*br - ci*bi;	
		buf[ix*2+1] = cr*bi + ci*br;	
		
		for(ix=0;ix<nx;ix++)
		{
			io[iw][ix*2]=buf[ix*2];
			io[iw][ix*2+1]=buf[ix*2+1];
		}
	}
	for(ix=0;ix<nx;ix++)
	{
		io[0][ix*2]=0.0;
		io[0][ix*2+1]=0.0;
	}
}

void sf_pscefd_step3(int iz,float **io) 
/*< step in depth >*/
{
	int ic;
	for(ic=0;ic<cas;ic++)
		pscefd_step3(iz,io);
}


void sf_pscefd_exit()
/*< free allocated storage >*/
{
	free(buf);
}

