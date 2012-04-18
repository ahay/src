/* 2-D Fourier finite-difference wave extrapolation, smooth point source, depress high frequency */
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
#include <math.h>
#include <limits.h>
#include "abcpass.h"
#include "ffdstep.h"
#include "srcsm.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nz, nt, ix, iz, it, nbt, nbb, nbl, nbr, nxb, nzb, isx, isz;
    float dt, dx, dz, o1, o2;
    float **old,  **cur,  **tmp, *wav;
    float  **v, v0, ***aa, w, g1, g2, ct, cb, cl, cr; /* top, bottom, left, right */
    float ax, az, factor;
    sf_file out, vel, source;
    bool opt;    /* optimal padding */
     
    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dz)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(vel,"o2",&o2)) o1=0.0;
    /*  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input"); */
    /*  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input"); */
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    if (!sf_getint("nbt",&nbt)) nbt=44;
    if (!sf_getint("nbb",&nbb)) nbb=44;
    if (!sf_getint("nbl",&nbl)) nbl=44;
    if (!sf_getint("nbr",&nbr)) nbr=44;

    if (!sf_getfloat("ct",&ct)) ct = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cl",&cl)) cl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cr",&cr)) cr = 0.01; /*decaying parameter*/

    if (!sf_getfloat("ax",&ax)) ax= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("az",&az)) az= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("factor",&factor)) factor= 5.0/6.0; /*suppress HF parameter*/

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
/*    sf_putfloat(out,"o1",x0); */
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o1",o1); 
    sf_putfloat(out,"o2",o2); 
    sf_putfloat(out,"o3",0.0); 

    nxb = nx + nbl + nbr;
    nzb = nz + nbt + nbb;
/*
    nkx = nxb;
    nkz = nzb;
    dkx = 1./(2.0*kiss_fft_next_fast_size(nxb-1)*dx);
    dkz = 1./(2.0*kiss_fft_next_fast_size(nzb-1)*dz);
*/

    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nxb,nzb);
    cur    =  sf_floatalloc2(nxb,nzb);
    aa     =  sf_floatalloc3(3,nxb,nzb);
    
    bd_init(nx,nz,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    v = sf_floatalloc2(nxb,nzb);
    /*input & extend velocity model*/
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(v[iz]+nbl,nx,vel);
         for (ix=0; ix<nbl; ix++){
             v[iz][ix] = v[iz][nbl];
         }
         for (ix=0; ix<nbr; ix++){
             v[iz][nx+nbl+ix] = v[iz][nx+nbl-1];
         }     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            v[iz][ix] = v[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            v[nz+nbt+iz][ix] = v[nz+nbt-1][ix];
        }
    }

    v0 =0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            v0 += v[iz][ix]*v[iz][ix];
         }
    }
    v0 = sqrtf(v0/(nxb*nzb));
 
    for (iz=0; iz < nzb; iz++){
         for (ix=0; ix < nxb; ix++) {
         w = v[iz][ix]*v[iz][ix];
         g1 = dt*dt*(v[iz][ix]*v[iz][ix]-v0*v0)/(12.0*dx*dx);
         g2 = dt*dt*(v[iz][ix]*v[iz][ix]-v0*v0)/(12.0*dz*dz);
         aa[iz][ix][1] = w*g1;
         aa[iz][ix][2] = w*g2;
         aa[iz][ix][0] = w-2.0*aa[iz][ix][1]-2.0*aa[iz][ix][2];
        }
      }

    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            cur[iz][ix] = 0.0;
            old[iz][ix] =  0.0; 
        }
    }

    /* propagation in time */
    ffdstep_init(nxb,nzb,dx,dz);
    srcsm_init(dz,dx);

    for (it=0; it < nt; it++) {

        ffdstep_dehf(old, cur, aa, nxb, nzb, v0, dt, ax, az, factor); 
        old[isz+nbt][isx+nbl] += wav[it];
        source_smooth(old,isz+nbt,isx+nbl,wav[it]);
        bd_decay(old); 
        bd_decay(cur); 
        tmp = old;
        old = cur;
        cur = tmp;
        
        for (iz=nbt; iz<nz+nbt; iz++){
             sf_floatwrite(cur[iz]+nbl,nx,out);
         }  
    }

    ffdstep_close();
    bd_close();
    free(**aa);
    free(*aa);
    free(aa);
    free(*v);     
    free(*cur);     
    free(*old);     
    free(v);     
    free(cur);     
    free(old);     
    exit(0); 
}           
           
