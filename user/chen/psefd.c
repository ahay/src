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

static int nw,nx;
static float df,dz,dx2;
static float **vel;
static sf_complex *buf;

void sf_psefd_init(int nx0,int nw0,
	float dx,float dz0,
	float f00,float **v0)
/*< initialize >*/
{
	nx=nx0;
	nw=nw0;

	dz=dz0;
	df=2.0*SF_PI*f00;
	vel=v0;

	dx2=(dx*dx);
	buf=(sf_complex*)sf_complexalloc(nx);
}

void sf_psefd_step3(int iz, sf_complex **io) 
/*< step in depth >*/
{
#ifdef SF_HAS_COMPLEX_H			// sf_complex = float complex
	int iw,ix;
	float t1, t2;
	sf_complex a0,a1,a2,b;
	//	0	a0i		1	a1i		0	a2i
//	a0r=0.0;
//	a1r=1.0;
//	a2r=0.0;

	for(iw=1; iw<nw; iw++)  // frequency slice
	{
		t1=dz/(iw*dx2*df);
		// ix=0 boundary
		ix=0;
		a1 = sf_cmplx(1.0,-t1*vel[ix][iz]);
		a2 = sf_cmplx(0,0.5*t1*vel[ix+1][iz]);
		t2 = df*dz*iw/vel[ix][iz];
		b = sf_cmplx(cos(t2),sin(t2));
		buf[ix] = b * (io[ix][iw]*a1 + io[ix+1][iw]*a2);

		for(ix=1;ix<nx-1;ix++)
		{
		    a0 = sf_cmplx(0.0,0.5*t1*vel[ix-1][iz]);
		    a1 = sf_cmplx(1.0,-t1*vel[ix][iz]);
		    a2 = sf_cmplx(0.0,0.5*t1*vel[ix+1][iz]);
			t2 = df*dz*iw/vel[ix][iz];
			b = sf_cmplx(cos(t2),sin(t2));
			buf[ix] = b * (io[ix-1][iw]*a0 + io[ix][iw]*a1 + io[ix+1][iw]*a2);
		}

		//ix=nx-1 boundary
		a0 = sf_cmplx(0.,0.5*t1*vel[ix-1][iz]);
		a1 = sf_cmplx(1.0,-t1*vel[ix][iz]);
		t2 = df*dz*iw/vel[ix][iz];
		b = sf_cmplx(cos(t2),sin(t2));
		buf[ix] = b * (io[ix-1][iw]*a0 + io[ix][iw]*a1);
		
		for(ix=0;ix<nx;ix++)
			io[ix][iw] = buf[ix];
	}
//	for(ix=0;ix<nx;ix++)
//		io[ix][0] = 0.0;
#else
	int iw,ix;
	float t1,a0i,a1i,a2i, br,bi,cr,ci;
	//	0	a0i		1	a1i		0	a2i
//	a0r=0.0;
//	a1r=1.0;
//	a2r=0.0;

	for(iw=1;iw<nw;iw++)  // frequency slice
	{
		t1=dz/(iw*dx2*df);
		// ix=0 boundary
		ix=0;
		a1i=-t1*vel[ix][iz];
		a2i=0.5*t1*vel[ix+1][iz];
		br=df*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[ix][iw].r
			- io[ix][iw].i   * a1i
			- io[ix+1][iw].i * a2i;
		ci=	+ io[ix][iw].i
			+ io[ix][iw].r   * a1i
			+ io[ix+1][iw].r * a2i;
		buf[ix].r = cr*br - ci*bi;	
		buf[ix].i = cr*bi + ci*br;		

		for(ix=1;ix<nx-1;ix++)
		{
			a0i=0.5*t1*vel[ix-1][iz];
			a1i=-t1*vel[ix][iz];
			a2i=0.5*t1*vel[ix+1][iz];
			br=df*dz*iw/vel[ix][iz];
			bi=sin(br);
			br=cos(br);
			cr=	  io[ix][iw].r 
				- io[ix-1][iw].i * a0i
				- io[ix][iw].i   * a1i
				- io[ix+1][iw].i * a2i;
			ci=	+ io[ix][iw].i
				+ io[ix-1][iw].r * a0i
				+ io[ix][iw].r   * a1i
				+ io[ix+1][iw].r * a2i;
			buf[ix].r = cr*br - ci*bi;	
			buf[ix].i = cr*bi + ci*br;		
		}

		//ix=nx-1 boundary
		a0i=0.5*t1*vel[ix-1][iz];
		a1i=-t1*vel[ix][iz];
		br=df*dz*iw/vel[ix][iz];
		bi=sin(br);
		br=cos(br);
		cr=	+ io[ix][iw].r
			- io[ix-1][iw].i * a0i
			- io[ix][iw].i   * a1i;
		ci=	+ io[ix][iw].i 
			+ io[ix-1][iw].r * a0i
			+ io[ix][iw].r   * a1i;
		buf[ix].r = cr*br - ci*bi;	
		buf[ix].i = cr*bi + ci*br;	
		
		for(ix=0;ix<nx;ix++)
		{
			io[ix][iw].r = buf[ix].r;
			io[ix][iw].i = buf[ix].i;
		}
	}
	for(ix=0;ix<nx;ix++)
	{
		io[ix][0].r = 0.0;
		io[ix][0].i = 0.0;
	}
#endif
}

void sf_psefd_close()
/*< release allocated memory >*/
{
	free(buf);
}

