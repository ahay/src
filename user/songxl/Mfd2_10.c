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
#include "srcsm.h"
int main(int argc, char* argv[]) 
{
    int nx, nt, ix, it, nz, iz, isx, isz;
    float dt, dx, dz, ox, oz, dx2, dz2, dt2, v2;
    float **nxt,  **old,  **cur, *wav, **v;
    float a, b1, b2, b3, b4, b5, c1, c2, c3, c4, c5; 
    sf_file out, vel, source;

    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"o1",&oz)) oz=0.0;
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o2",&ox)) ox=0.0;
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    sf_putint(out,"n1",nz);
    sf_putfloat(out,"d1",dz);
    sf_putfloat(out,"o1",oz); 
    sf_putint(out,"n2",nx);
    sf_putfloat(out,"d2",dx);
    sf_putfloat(out,"o2",ox); 
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
    v      =  sf_floatalloc2(nz,nx);

    sf_floatread(v[0],nz*nx,vel);

    dx2 = dx*dx;
    dz2 = dz*dz;
    dt2 = dt*dt;
    b1  = 5.0/(3.0*dz2)*dt2;
    b2  = -5.0/(21.0*dz2)*dt2;
    b3  = 5.0/(126.0*dz2)*dt2;
    b4  = -5.0/(1008.0*dz2)*dt2;
    b5  = 1.0/(3150.0*dz2)*dt2;
    c1  = 5.0/(3.0*dx2)*dt2;
    c2  = -5.0/(21.0*dx2)*dt2;
    c3  = 5.0/(126.0*dx2)*dt2;
    c4  = -5.0/(1008.0*dx2)*dt2;
    c5  = 1.0/(3150.0*dx2)*dt2;
    a   = -2.0*(b1+b2+b3+b4+b5+c1+c2+c3+c4+c5);

    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            cur[ix][iz] = 0.0;
            old[ix][iz] = 0.0; 
        }
    }
    cur[isx][isz] = wav[0];
    srcsm_init(dx,dz);
    source_smooth(cur,isx,isz,wav[0]);

    /* propagation in time */
    for (it=0; it < nt; it++) {

        sf_floatwrite(cur[0],nz*nx,out);
        for (ix=0; ix < nx; ix++) {
            for (iz=0; iz < nz; iz++) {
                  nxt[ix][iz] = 0.0; 
                }
         }  




	 for (ix=5; ix < nx-5; ix++) {  
	     for (iz=5; iz < nz-5; iz++) {  
                 v2 = v[ix][iz]*v[ix][iz];
                 nxt[ix][iz]  = cur[ix][iz]*a*v2
                              + (cur[ix][iz-1]+cur[ix][iz+1])*b1*v2
                              + (cur[ix][iz-2]+cur[ix][iz+2])*b2*v2
                              + (cur[ix][iz-3]+cur[ix][iz+3])*b3*v2
                              + (cur[ix][iz-4]+cur[ix][iz+4])*b4*v2
                              + (cur[ix][iz-5]+cur[ix][iz+5])*b5*v2
                              + (cur[ix-1][iz]+cur[ix+1][iz])*c1*v2
                              + (cur[ix-2][iz]+cur[ix+2][iz])*c2*v2
                              + (cur[ix-3][iz]+cur[ix+3][iz])*c3*v2
                              + (cur[ix-4][iz]+cur[ix+4][iz])*c4*v2
                              + (cur[ix-5][iz]+cur[ix+5][iz])*c5*v2
                              - old[ix][iz] + 2.0*cur[ix][iz];

/*
                 nxt[ix][iz]  = 
                               0.5*(cur[ix][iz-1]+cur[ix][iz+1]-2.0*cur[ix][iz])*b1[ix][iz]
                              + 0.5*(cur[ix][iz-2]+cur[ix][iz+2]-2.0*cur[ix][iz])*b2[ix][iz]
                              + 0.5*(cur[ix-1][iz]+cur[ix+1][iz]-2.0*cur[ix][iz])*c1[ix][iz]
                              + 0.5*(cur[ix-2][iz]+cur[ix+2][iz]-2.0*cur[ix][iz])*c2[ix][iz]
                              + 0.5*(cur[ix+1][iz+1]+cur[ix-1][iz-1]-2.0*cur[ix][iz])*d1[ix][iz]
                              + 0.5*(cur[ix+1][iz-1]+cur[ix-1][iz+1]-2.0*cur[ix][iz])*d2[ix][iz]
                              + 2.0*cur[ix][iz] - old[ix][iz];
*/
             }
         }  
	 for (ix=5; ix < nx-5; ix++) {  
	     for (iz=5; iz < nz-5; iz++) {  
                 old[ix][iz] = cur[ix][iz];
                 cur[ix][iz] = nxt[ix][iz];
             }
         }  
         cur[isx][isz] += wav[it];
         source_smooth(cur,isx,isz,wav[it]);
         
                 
    }
    exit(0); 
}           
           
