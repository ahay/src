/* 2-D image ray tracing using HWT */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
    int nt, nx, nz, it, ix, iz, order;
    float dt, dx, dz;
    float **vv, *xx, *xp, *zz, *zp, *vd;
    sf_eno2 vmap;
    sf_file vel, dix;

    sf_init(argc,argv);

    vel = sf_input("in");
    dix = sf_output("out");

    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");

    if (!sf_histint(vel,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dz)) sf_error("No d2= in input");
 
    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    
    sf_putint(dix,"n2",3);

    sf_putint(dix,"n3",nt);
    sf_putfloat(dix,"d3",dt);
    sf_putfloat(dix,"o3",0.);
    sf_putstring(dix,"label3","Time");
    sf_putstring(dix,"unit3","s");

    vv = sf_floatalloc2(nx,nz);
    sf_floatread(vv[0],nz*nx,vel);

    if (!sf_getint("order",&order)) order=3;
    /* interpolation order */

    vmap = sf_eno2_init(order,nx,nz);
    sf_eno2_set(vmap,vv);

    xx = sf_floatalloc(nx);
    xp = sf_floatalloc(nx);

    zz = sf_floatalloc(nx);
    zp = sf_floatalloc(nx);
  
    vd = sf_floatalloc(nx);

    exit(0);
}
