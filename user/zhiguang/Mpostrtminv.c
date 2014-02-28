/* 2-D least-squares RTM of incomplete zero-offset data */
/*
 Copyright (C) 2014 University of Texas at Austin
 
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
#include "postrtm.h"

int main(int argc, char* argv[])
{
    bool verb;
    int ix, iz;
    int nx, nz, nt, n0, nxz, nxt, padx, padz, n2, padnx, padnz, niter;
    float dx, dz, dt, dt2;
    
    int *head;
    float *mm, *dd;
    float **vv, **padvv;
    
    sf_file in, out, vel, mask;
    sf_init(argc, argv);
    
    in=sf_input("in");
    out=sf_output("out");
    vel=sf_input("velocity");
    mask=sf_input("mask");
    
    /* Dimensions */
    if(!sf_histint(vel, "n1", &nz)) sf_error("No n1= in velocity");
    if(!sf_histint(vel, "n2", &nx)) sf_error("No n2= in velocity");
    if(!sf_histfloat(vel, "d1", &dz)) sf_error("No d1= in velocity");
    if(!sf_histfloat(vel, "d2", &dx)) sf_error("No d2= in velocity");
    
    if(!sf_histint(in, "n1", &nt)) sf_error("No n1= in data");
    if(!sf_histfloat(in, "d1", &dt)) sf_error("No d1= in data");
    if(!sf_histint(in, "n2", &n2) || n2!=nx) sf_error("Need n2=%d in data", nx);
    if(!sf_histint(mask, "n1", &n2) || n2!=nx) sf_error("Need n1=%d in mask", nx);
    
    sf_putint(out, "n1", nz);
    sf_putfloat(out, "d1", dz);
    sf_putfloat(out, "o1", 0.0);
    sf_putstring(out, "label1", "Depth");
    sf_putstring(out, "unit1", "km");
    sf_putstring(out, "label2", "Lateral");
    sf_putstring(out, "unit2", "km");
    
    if(!sf_getbool("verb", &verb)) verb=true;
    if(!sf_getint("niter", &niter)) niter=5;
    if(!sf_getint("padx", &padx)) padx=nz/2;
    if(!sf_getint("padz", &padz)) padz=nz/2;
    if(!sf_getint("n0", &n0)) n0=0;
    /* surface */
    
    padnx=nx+2*padx;
    padnz=nz+2*padz;
    n0=padz+n0;
    nxz=nx*nz;
    nxt=nx*nt;
    
    /* allocate arrays */
    vv=sf_floatalloc2(nz, nx);
    padvv=sf_floatalloc2(padnz, padnx);
    dd=sf_floatalloc(nxt);
    mm=sf_floatalloc(nxz);
    head=sf_intalloc(nx);
    
    /* read velocity, mask and data */
    sf_floatread(vv[0], nxz, vel);
    sf_intread(head, nx, mask);
    sf_floatread(dd, nxt, in);
    
    /* pad boundary */
    dt2=dt*dt;
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            padvv[ix+padx][iz+padz]=vv[ix][iz]*vv[ix][iz]*dt2;
    for(iz=0; iz<padz; iz++){
        for(ix=padx; ix<nx+padx; ix++){
            padvv[ix][iz]=padvv[ix][padz];
            padvv[ix][iz+nz+padz]=padvv[ix][nz+padz-1];
        }
    }
    for(ix=0; ix<padx; ix++){
        for(iz=0; iz<padnz; iz++){
            padvv[ix][iz]=padvv[padx][iz];
            padvv[ix+nx+padx][iz]=padvv[nx+padx-1][iz];
        }
    }
    
    postrtm_init(nx, nz, nt, n0, padx, padz, padnx,
                    padnz, dx, dz, head, padvv);
    
    sf_solver(postrtm_lop, sf_cgstep, nxz, nxt, mm, dd, niter, "verb", verb, "end");
    
    sf_floatwrite(mm, nxz, out);
    
    exit(0);
}
