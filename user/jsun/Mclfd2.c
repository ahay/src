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

int main(int argc, char* argv[]) 
{
    int nx, nt, ix, it, nz, iz, isx, isz, nxz, ik;
    int LEN,SIZE;
    int *stmp, *s1, *s2;
    float dt, dx, dz;
    sf_complex **nxt, **old, **cur, *src;
    sf_complex ***coef, plap;
    sf_file out, vel, source, G;

    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    source = sf_input("in");   /* source wavlet*/
    G = sf_input("G");   /* source wavlet*/

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_COMPLEX != sf_gettype(source)) sf_error("Need complex source");
    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");
    if (!sf_getint("size",&SIZE)) sf_error("Need size");
    if (!sf_histint(G,"n1",&nxz)) sf_error("No n1= in input");
    if (!sf_histint(G,"n2",&LEN)) sf_error("No n2= in input");
    if (nx*nz != nxz) sf_error("nx*nz != nxz");

    sf_putint(out,"n1",nz);
    sf_putfloat(out,"d1",dz);
    sf_putint(out,"n2",nx);
    sf_putfloat(out,"d2",dx);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o3",0.0); 
    sf_settype(out,SF_COMPLEX);

    sf_warning("nt=%d",nt);
    sf_warning("dt=%g",dt);
    src =  sf_complexalloc(nt);
    sf_complexread(src,nt,source);

    old  = sf_complexalloc2(nz,nx);
    cur  = sf_complexalloc2(nz,nx);
    nxt  = sf_complexalloc2(nz,nx);

    coef = sf_complexalloc3(nz,nx,LEN);
    
    sf_complexread(coef[0][0],nz*nx*LEN,G);

    stmp = sf_intalloc(SIZE);
    s1 = sf_intalloc(LEN);
    s2 = sf_intalloc(LEN);

    for (int ix=0; ix<SIZE; ix++) stmp[ix]= ix - (SIZE-1)/2;

    ik = 0;
    for (int ix=0; ix<SIZE; ix++){
        for (int iz=0; iz<SIZE; iz++){
            if((stmp[ix] == 0) || (stmp[iz] == 0)) {
                s1[ik]=stmp[iz];
                s2[ik]=stmp[ix];
                ik++;
            }
        }
    }
    if (ik!=LEN) sf_error("Stencil size mismatch!!!");
    
    for (ix=0; ix < nx; ix++) {
        for (iz=0; iz < nz; iz++) {
            cur[ix][iz] = sf_cmplx(0,0);
            old[ix][iz] = sf_cmplx(0,0);
            nxt[ix][iz] = sf_cmplx(0,0);
        }
    }

    /* propagation in time */
    for (it=0; it < nt; it++) {
        sf_warning("it=%d;",it);

#ifdef SF_HAS_COMPLEX_H
        cur[isx][isz] += src[it];
#else
        cur[isx][isz] = sf_cadd(cur[isx][isz],src[it]);
#endif

        sf_complexwrite(cur[0],nz*nx,out);

        for (ix=2; ix < nx-2; ix++) {  
            for (iz=2; iz < nz-2; iz++) {  

                // apply lowrank FD coefficients
                plap = sf_cmplx(0,0);
                for (ik=0; ik < LEN; ik++) {
#ifdef SF_HAS_COMPLEX_H
                    plap += coef[ik][ix][iz]*cur[ix+s2[ik]][iz+s1[ik]];
#else
                    plap = sf_cadd(plap,coef[ik][ix][iz]*cur[ix+s2[ik]][iz+s1[ik]]);
#endif
                }

                nxt[ix][iz] = plap; 
                
            }
        }  

        for (ix=2; ix < nx-2; ix++) {  
            for (iz=2; iz < nz-2; iz++) {  
                old[ix][iz] = cur[ix][iz];
                cur[ix][iz] = nxt[ix][iz];
            }
        }  

    }
    sf_warning(".");

    exit(0); 
}           
           
