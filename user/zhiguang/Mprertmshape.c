/* 2D prestack FD RTM with LS and shaping regularization */
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
#include "preshape.h"

int main(int argc, char* argv[])
{
    bool verb;
    int ix, iz, nd, nm;
    int nz, nx, nt, nr, ns, nw, nsource, dsource, ndelay;
    int padx, padz, padnx, padnz, niter, rect1, rect2, order;
    int dr_v, ds_v, r0_v, s0_v, zr_v, zs_v; 
    float dz, dx, dt, dr, ds, z0, x0, t0, r0, s0;
    float zr, zs, padx0, padz0, tdelay, lambda, dt2;
    
    float **vv, **padvv, **dp, *ww, *mm, *dd;
    
    sf_file in, out, vel, wavelet, dip;
    sf_init(argc, argv);
    
    in=sf_input("in");
    out=sf_output("out");
    vel=sf_input("velocity");
    wavelet=sf_input("wavelet");
    dip=sf_input("dip");
    
    /* Dimensions */
    if(!sf_histint(vel, "n1", &nz)) sf_error("No n1= in velocity");
    if(!sf_histfloat(vel, "d1", &dz)) sf_error("No d1= in velocity");
    if(!sf_histfloat(vel, "o1", &z0)) sf_error("No o1= in velocity");
    
    if(!sf_histint(vel, "n2", &nx)) sf_error("No n2= in velocity");
    if(!sf_histfloat(vel, "d2", &dx)) sf_error("No d2= in velocity");
    if(!sf_histfloat(vel, "o2", &x0)) sf_error("No o2= in velocity");
    
    if(!sf_histint(wavelet, "n1", &nt)) sf_error("No n1= in wavelet");
    if(!sf_histfloat(wavelet, "d1", &dt)) sf_error("No d1= in wavelet");
    if(!sf_histfloat(wavelet, "o1", &t0)) sf_error("No o1= in wavelet");
    
    if(!sf_histint(in, "n2", &nr)) sf_error("No n2= in input");
    if(!sf_histfloat(in, "d2", &dr)) sf_error("No d2= in input");
    if(!sf_histfloat(in, "o2", &r0)) sf_error("No o2= in input");
        
    if(!sf_histint(in, "n3", &ns)) sf_error("No n3= in input");
    if(!sf_histfloat(in, "d3", &ds)) sf_error("No d3= in input");
    if(!sf_histfloat(in, "o3", &s0)) sf_error("No o3= in input");
        
    sf_putint(out, "n1", nz);
    sf_putfloat(out, "d1", dz);
    sf_putfloat(out, "o1", z0);
    sf_putint(out, "n2", nx);
    sf_putfloat(out, "d2", dx);
    sf_putfloat(out, "o2", x0);
    sf_putint(out, "n3", 1);
    sf_putstring(out, "label1", "Depth");
    sf_putstring(out, "unit1", "km");
    sf_putstring(out, "label2", "Lateral");
    sf_putstring(out, "unit2", "km");
    
    if(!sf_getbool("verb", &verb)) verb=true;
    if(!sf_getfloat("zr", &zr)) zr=0.0;
    if(!sf_getfloat("zs", &zs)) zs=0.0;
    if(!sf_getint("padz", &padz)) padz=nz;
    if(!sf_getint("padx", &padx)) padx=nz;
    if(!sf_getint("nw", &nw)) nw=1;
    if(!sf_getint("nsource", &nsource)) nsource=1;
    if(!sf_getint("dsource", &dsource)) dsource=0;
    if(!sf_getfloat("tdelay", &tdelay)) tdelay=0;
    
    if(!sf_getint("rect1", &rect1)) rect1=3;
    if(!sf_getint("rect2", &rect2)) rect2=3;
    /* smoothing radius */
    if(!sf_getint("order", &order)) order=1;
    /* accuracy order */
    if(!sf_getfloat("lambda", &lambda)) lambda=1;
    /* operator scaling for inversion */
    if(!sf_getint("niter", &niter)) niter=3;
    
    padnx=nx+2*padx;
    padnz=nz+2*padz;
    padx0=x0-dx*padx;
    padz0=z0-dz*padz;
    
    dr_v=(dr/dx)+0.5;
    r0_v=(r0-padx0)/dx+0.5;
    zr_v=(zr-padz0)/dz+0.5;
    
    ds_v=(ds/dx)+0.5;
    s0_v=(s0-padx0)/dx+0.5;
    zs_v=(zs-padz0)/dz+0.5;
    
    nd=nt*nr*ns;
    nm=nz*nx;
    ndelay=tdelay/dt;
    
    dd=sf_floatalloc(nd);
    mm=sf_floatalloc(nm);
    dp=sf_floatalloc2(nz, nx);
    vv=sf_floatalloc2(nz, nx);
    padvv=sf_floatalloc2(padnz, padnx);
    ww=sf_floatalloc(nt);
    
    sf_floatread(ww, nt, wavelet);
    sf_floatread(vv[0], nm, vel);
    sf_floatread(dd, nd, in);
    sf_floatread(dp[0], nm, dip);
    
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

    sf_warning("nx=%d nz=%d nt=%d nr=%d ns=%d", nx,  nz, nt, nr, ns);
    sf_warning("padx=%d padnx=%d padz=%d padnz=%d padx0=%.3f padz0=%.3f", padx, padnx, padz, padnz, padx0, padz0);
    sf_warning("dx=%.3f dz=%.3f dt=%.3f dr=%.3f ds=%.3f", dx, dz, dt, dr, ds);
    sf_warning("x0=%.3f z0=%.3f t0=%.3f r0=%.3f s0=%.3f", x0, z0, t0, r0, s0);
    sf_warning("dr_v=%d r0_v=%d zr_v=%d ds_v=%d s0_v=%d zs_v=%d", dr_v, r0_v, zr_v, ds_v, s0_v, zs_v);
    sf_warning("nsource=%d dsource=%d tdelay=%.3f ndelay=%d", nsource, dsource, tdelay, ndelay);
    
    preshape_init(verb, nz, nx, nt, nr, ns, nw, nsource, dsource, ndelay,
                  dx, dz, padx, padz, padnx, padnz,
                  nm, nd, dr_v, ds_v, r0_v, s0_v, zr_v, zs_v,
                  padvv, ww, rect1, rect2, order, dp, lambda);
                                      
    preshape(niter, mm, dd);
    preshape_close();
        
    sf_floatwrite(mm, nm, out);
    
    exit(0);
}
