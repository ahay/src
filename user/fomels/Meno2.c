/* Convert velocity to slowness and compute gradient using ENO interpolation */
/*
  Copyright (C) 2009 University of Texas at Austin

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

int main(int argc, char* argv[])
{
    int nz, nx, iz, ix, order;
    float dz, dx, s;
    float **v, **vx, **vz, grad[2];
    bool is_inverse;
    sf_eno2 slow;
    sf_file vel, out;

    sf_init(argc,argv);
    vel = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");

    if (!sf_histint(vel,"n1",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d2= in input");

    sf_putint(out,"n3",3);

    if (!sf_getint("order",&order)) order=3;
    /* interpolation order */

    if (!sf_getbool("is_inverse",&is_inverse)) is_inverse=1;
    /* make vel to slowness */


    slow = sf_eno2_init (order,nz,nx);
    v = sf_floatalloc2(nz,nx);

    sf_floatread(v[0],nz*nx,vel);

    /* convert velocity to slowness */
    if (is_inverse) {
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		v[ix][iz] = 1./v[ix][iz];
	    }
	}
    }
    sf_eno2_set(slow,v);

    vz = sf_floatalloc2(nz,nx);
    vx = sf_floatalloc2(nz,nx);

    for (ix=0; ix < nx; ix++) { 
	for (iz=0; iz < nz; iz++) {
	    sf_eno2_apply (slow,iz,ix,0.,0.,&s,grad,BOTH);
	    v[ix][iz] = s;
	    vz[ix][iz] = grad[0]/dz;
	    vx[ix][iz] = grad[1]/dx;
	}
    }

    sf_floatwrite(v[0],nz*nx,out);
    sf_floatwrite(vz[0],nz*nx,out);
    sf_floatwrite(vx[0],nz*nx,out);
    
    exit(0);
}
