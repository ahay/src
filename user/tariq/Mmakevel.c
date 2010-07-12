/* Make a velocity function v(x,y,z)*/
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

#include <math.h>

#include <rsf.h>

#include "frannor.h"

#define PI 3.14159265359



/*
 * Extracted from CWP Author: Dave Hale
 *
 */

int
main (int argc, char **argv)
{
	int nx,ny,nz,ix,iy,iz;
	float dx,dy,dz,fx,fy,fz,x,y,z,v000,dvdx,dvdy,dvdz,
		xlens,ylens,zlens,dlens,tlens,vlens,xn,ynn,zn,
		/* abot,bbot,atop,btop, */
		vran,vzran, 
		vzc,z1c,z2c,l1c,l2c,exc,ac,bc,vtemp,
		*v,*vz;
	/* float zbot,ztop;*/
	char *vzfile;
	sf_file out, vzf;

	/* hook up getpar to handle the parameters */
	sf_init (argc, argv);

	out = sf_output("out");

	if (NULL != (vzfile = sf_getstring("vzfile"))) {
	  vzf = sf_input(vzfile);
	  free(vzfile);
	} else {
	  vzf = NULL;
	}
	
	/* get required parameters */
	if (!sf_getint("n2",&nx)) sf_error("must specify n2!\n");
	/* number of x samples (2nd dimension), must be provided!*/
	if (!sf_getint("n1",&nz)) sf_error("must specify n1!\n");
	/* number of z samples (1st dimension)), must be provided!*/
	
	/* get optional parameters */
	if (!sf_getint("n3",&ny)) ny = 1;
	/* number of y samples (3rd dimension)*/
	if (!sf_getfloat("d2",&dx)) dx = 1.0;
	/* 2nd dimension sampling interval*/
	if (!sf_getfloat("d3",&dy)) dy = 1.0;
	/* 3rd dimension sampling interval*/
	if (!sf_getfloat("d1",&dz)) dz = 1.0;
	/* 1st dimension sampling interval*/
	if (!sf_getfloat("o2",&fx)) fx = 0.0;
	/* Origin 2nd dimension*/
	if (!sf_getfloat("o3",&fy)) fy = 0.0;
	/* Origin 3rd dimension*/
	if (!sf_getfloat("o1",&fz)) fz = 0.0;
	/* Origin 1st dimension*/
	if (!sf_getfloat("v000",&v000)) v000 = 2.0;
	/* velocity at (x=0,y=0,z=0)*/
	if (!sf_getfloat("dvdx2",&dvdx)) dvdx = 0.0;
	/* velocity gradient with respect to 2nd dimension*/
	if (!sf_getfloat("dvdx3",&dvdy)) dvdy = 0.0;
	/* velocity gradient with respect to 3rd dimension*/
	if (!sf_getfloat("dvdx1",&dvdz)) dvdz = 0.0;
	/* velocity gradient with respect to 1st dimension*/
	if (!sf_getfloat("x2lens",&xlens)) xlens = fx;
	/* 2nd dimension coordinate of center of parabolic lens*/
	if (!sf_getfloat("x3lens",&ylens)) ylens = fy;
	/* 3rd dimension coordinate of center of parabolic lens*/
	if (!sf_getfloat("x1lens",&zlens)) zlens = fz;
	/* 1st dimension coordinate of center of parabolic lens*/
	if (!sf_getfloat("vlens",&vlens)) vlens = 0.0;
	/* velocity perturbation in parabolic lens*/
	if (!sf_getfloat("dlens",&dlens)) dlens = 1.0;
	/* diameter of parabolic lens*/
	if (!sf_getfloat("tlens",&tlens)) tlens = 1.0;
	/* thickness of parabolic lens*/
	if (!sf_getfloat("vran",&vran)) vran = 0.0;
	/* standard deviation of random perturbation*/
	vzfile = sf_getstring("vx1file");
	/* file containing v(z) 1st dimension profile*/
	if (!sf_getfloat("vx1ran",&vzran)) vzran = 0.0;
	/* standard deviation of random perturbation to 1st dimension*/
	if (!sf_getfloat("vx1c",&vzc)) vzc = 0.0;
	/* 1st dimension v(z) chirp amplitude*/
	if (!sf_getfloat("x11c",&z1c)) z1c = fz;
	/* 1st dimension at which to begin chirp*/
	if (!sf_getfloat("x12c",&z2c)) z2c = fz+(nz-1)*dz;
	/* 1st dimension at which to end chirp*/
	if (!sf_getfloat("l1c",&l1c)) l1c = dz;
	/* wavelength at beginning of chirp*/
	if (!sf_getfloat("l2c",&l2c)) l2c = dz;
	/*wavelength at end of chirp*/
	if (!sf_getfloat("exc",&exc)) exc = 1.0;
	/* exponent of chirp*/

	sf_setformat(out,"native_float");
	sf_putint(out,"n1",nz);
	sf_putint(out,"n2",nx);
	sf_putint(out,"n3",ny);
	sf_putfloat(out,"d1",dz);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"d3",dy);
	sf_putfloat(out,"o1",fz);
	sf_putfloat(out,"o2",fx);
	sf_putfloat(out,"o3",fy);
	
	/* compute chirp constants */
	bc = PI/(z2c-z1c)*(1.0/l2c-1.0/l1c);
	ac = 2.0*PI/l1c - 2.0*bc*z1c;
	
	/* allocate space */
	v = sf_floatalloc(nz);
	vz = sf_floatalloc(nz);
	
	if (NULL != vzf) {
	  sf_floatread(vz,nz,vzf);
	  sf_fileclose(vzf);
	} else {
	  for (iz=0; iz<nz; ++iz)
	    vz[iz] = 0.0;
	}

	/* random v(z) perturbation */
	/*for (iz=0; iz<nz; ++iz)
	  vz[iz] += vzran*frannor();*/


	/* loop over x */
	for (ix=0,x=fx; ix<nx; ++ix,x+=dx) {
	
		/* loop over y */
		for (iy=0,y=fy; iy<ny; ++iy,y+=dy) {
		
			/* compute top and bottom of lens */
/*
			ztop = atop+btop*(pow(x-xlens,2)+pow(y-ylens,2));
			zbot = abot+bbot*(pow(x-xlens,2)+pow(y-ylens,2));
*/
			
			/* loop over z */
			for (iz=0,z=fz; iz<nz; ++iz,z+=dz) {
				
				/* v(z) profile */
				v[iz] = vz[iz];
				
				/* constant + constant gradient */
				v[iz] += v000+x*dvdx+y*dvdy+z*dvdz;
				
				/* lens */
				xn = 2.0*(x-xlens)/dlens;
				ynn = 2.0*(y-ylens)/dlens;
				zn = 2.0*(z-zlens)/tlens;
				v[iz] += vlens*exp(-xn*xn-ynn*ynn-zn*zn);
				/*
				if (z>zbot && z<ztop) v[iz] += vlens;
				*/

				/* chirp */
				if (z>z1c && z<z2c) {
					vtemp = sin((ac+bc*z)*z);
					if (vtemp<0.0)
						v[iz] -= vzc*pow(-vtemp,exc);
					else
						v[iz] += vzc*pow(vtemp,exc);
				}

				/* random perturbation */
				/*v[iz] += vran*frannor();*/
				/*warn("v00=%f dvdx=%f dvdz=%f vz=%f",v000,dvdx,dvdz,v[iz]);*/
			}
			
			/* write velocity function */
			sf_floatwrite (v,nz,out);
		}
	}
	
	exit (0);
}
