/* 2-D Fourth-order Optimized Finite-difference wave extrapolation */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <math.h>
#include <limits.h>
#include "abcpass.h"
int main(int argc, char* argv[]) 
{
    int nx, nt, ix, it, nz, iz, isx, isz, nxz, na;
    float dt, dx, dz;
    float **nxt,  **old,  **cur, *wav;
    float **a, **b1, **b2, **c1, **c2; 
    sf_file out, vel, source, G;

    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    source = sf_input("in");   /* source wavlet*/
    G = sf_input("G");   /* source wavlet*/

//    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");
    if (!sf_histint(G,"n1",&nxz)) sf_error("No n1= in input");
    if (!sf_histint(G,"n2",&na)) sf_error("No n2= in input");
    if (nx*nz != nxz) sf_error("nx*nz != nxz");

    sf_putint(out,"n1",nz);
    sf_putfloat(out,"d1",dz);
    sf_putint(out,"n2",nx);
    sf_putfloat(out,"d2",dx);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o3",0.0); 

    sf_warning("nt=%d",nt);
    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);
    sf_warning("dt=%g",dt);

    old    =  sf_floatalloc2(nz,nx);
    cur    =  sf_floatalloc2(nz,nx);
    nxt    =  sf_floatalloc2(nz,nx);

    a     =  sf_floatalloc2(nz,nx);
    b1    =  sf_floatalloc2(nz,nx);
    b2    =  sf_floatalloc2(nz,nx);
    c1    =  sf_floatalloc2(nz,nx);
    c2    =  sf_floatalloc2(nz,nx);
    
    sf_floatread(a[0],nz*nx,G);
    sf_floatread(b1[0],nz*nx,G);
    sf_floatread(b2[0],nz*nx,G);
    sf_floatread(c1[0],nz*nx,G);
    sf_floatread(c2[0],nz*nx,G);

    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            cur[ix][iz] = 0.0;
            old[ix][iz] = 0.0; 
        }
    }
    cur[isx][isz] = wav[0];

    /* propagation in time */
    for (it=0; it < nt; it++) {

        sf_floatwrite(cur[0],nz*nx,out);
        for (ix=0; ix < nx; ix++) {
            for (iz=0; iz < nz; iz++) {
                  nxt[ix][iz] = 0.0; 
                }
         }  




	 for (ix=2; ix < nx-2; ix++) {  
	     for (iz=2; iz < nz-2; iz++) {  
                 nxt[ix][iz]  = cur[ix][iz]*(2.0*a[ix][iz])
                              + (cur[ix][iz-1]+cur[ix][iz+1])*b1[ix][iz]
                              + (cur[ix][iz-2]+cur[ix][iz+2])*b2[ix][iz]
                              + (cur[ix-1][iz]+cur[ix+1][iz])*c1[ix][iz]
                              + (cur[ix-2][iz]+cur[ix+2][iz])*c2[ix][iz]
                              + 2.0*cur[ix][iz] - old[ix][iz];
                 old[ix][iz] = cur[ix][iz];
                 cur[ix][iz] = nxt[ix][iz];
                 cur[isx][isz] += wav[it];
             }
         }  
         
                 
    }
    exit(0); 
}           
           
