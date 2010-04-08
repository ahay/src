/* Linearized complex eikonal equation */
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
    int dim, n[3], ix, iy, iz;
    float v0, b, o[3], d[3], l[3], source[3], s, xi, *t;
    sf_file vel, time;

    sf_init(argc,argv);
    vel = sf_input("in");
    time = sf_output("out");

    dim = sf_filedims(vel,n);

    if(!sf_histfloat(vel,"o1",&o[0])) o[0]=0.;
    if(!sf_histfloat(vel,"o2",&o[1])) o[1]=0.;
    if(!sf_histfloat(vel,"o3",&o[2])) o[2]=0.;

    if(!sf_histfloat(vel,"d1",&d[0])) sf_error("No d1 in input");
    if(!sf_histfloat(vel,"d2",&d[1])) sf_error("No d2 in input");
    if(!sf_histfloat(vel,"d3",&d[2])) d[2]=d[1];
    
    if(!sf_getfloat("s1",&source[0])) source[0]=o[0];
    if(!sf_getfloat("s2",&source[1])) source[1]=o[1];
    if(!sf_getfloat("s3",&source[2])) source[2]=o[2];

    if(!sf_getfloat("v0",&v0)) v0=1.;
    if(!sf_getfloat("b",&b)) b=0.;
    if(!sf_getfloat("s",&s)) s=0.;

    t = sf_floatalloc(n[2]*n[1]*n[0]);

    if (dim < 3)
	n[2] = 1;

    for (ix=0; ix < n[2]; ix++) {
	l[2] = ix*d[2]+o[2]-source[2];

	for (iy=0; iy < n[1]; iy++) {
	    l[1] = iy*d[1]+o[1]-source[1];

	    for (iz=0; iz < n[0]; iz++) {
		l[0] = iz*d[0]+o[0]-source[0];
		
		xi = b*b*(l[0]*l[0]+l[1]*l[1]+l[2]*l[2])/(2*(v0+b*l[0])*v0);
		t[iz+iy*n[0]+ix*n[1]*n[0]] = (1/b)*logf(1+xi+sqrtf((2+xi)*xi));
	    }
	}
    }

    sf_floatwrite(t,n[2]*n[1]*n[0],time);

    exit(0);
}
