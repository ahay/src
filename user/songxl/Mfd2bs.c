/* 2-D Fourth-order Finite-difference wave extrapolation */
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
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nt, ix, it, nz, iz, nb, nxb, nzb, isx, isz;
    float dt, dx, dz;
    float **nxt,  **old,  **cur,  **dercur, **derold, *wav;
    float  **v, **aa, *wb, c; 
    sf_file out, vel, source;
    bool opt;    /* optimal padding */
    /* #ifdef _OPENMP
   int nth;
   #endif */
     

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
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    if (!sf_getint("nb",&nb)) nb=30;

    if (!sf_getfloat("c",&c)) c = 0.01; /*decaying parameter*/
    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
/*    sf_putfloat(out,"o1",x0); */
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);
    sf_putint(out,"n3",nt);
    sf_putfloat(out,"d3",dt);
    sf_putfloat(out,"o3",0.0); 

    nxb = nx + 2*nb + 4;
    nzb = nz + 2*nb + 4;

    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nxb,nzb);
    cur    =  sf_floatalloc2(nxb,nzb);
    nxt    =  sf_floatalloc2(nxb,nzb);
    derold    =  sf_floatalloc2(nxb,nzb);
    dercur    =  sf_floatalloc2(nxb,nzb);
    aa     =  sf_floatalloc2(3,2);
    
    v = sf_floatalloc2(nxb-4,nzb-4);
    /*input & extend velocity model*/
    for (iz=nb; iz<nz+nb; iz++){
        sf_floatread(v[iz]+nb,nx,vel);
         for (ix=0; ix<nb; ix++){
             v[iz][ix] = v[iz][nb];
             v[iz][nx+nb+ix] = v[iz][nx+nb-1];
         }
    }     
    for (iz=0; iz<nb; iz++){
        for (ix=0; ix<nxb; ix++){
            v[iz][ix] = v[nb][ix];
            v[nz+nb+iz][ix] = v[nz+nb-1][ix];
        }
    }

 
    aa[0][0] = -30.0/(12.0*dx*dx); 
    aa[0][1] =  16.0/(12.0*dx*dx); 
    aa[0][2] = - 1.0/(12.0*dx*dx); 
    aa[1][0] = -30.0/(12.0*dz*dz); 
    aa[1][1] =  16.0/(12.0*dz*dz); 
    aa[1][2] = - 1.0/(12.0*dz*dz); 

    wb =  nb? sf_floatalloc(nb): NULL;
    abc_cal(0,nb,c,wb);


    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            cur[iz][ix] = 0.0;
        }
    }
    cur[isz+nb+2][isx+nb+2] = wav[0];
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            old[iz][ix] =  0.0; 
            derold[iz][ix] =cur[iz][ix]/dt;
           }
         }
    for (iz=nb+2; iz<nz+nb+2; iz++){
        sf_floatwrite(cur[iz]+nb+2,nx,out);
    }
/*
    #ifdef _OPENMP
    #pragma omp parallel
   {nth = omp_get_num_threads();}
    sf_warning("using %d threads",nth);
    #endif
*/
    /* propagation in time */
    for (it=1; it < nt; it++) {

         for (iz=0; iz < nzb; iz++){
             for (ix=0; ix < nxb; ix++){ 
                  nxt[iz][ix] = 0.0; 
                }
         }  



	 for (iz=2; iz < nzb-2; iz++) {  
	     for (ix=2; ix < nxb-2; ix++) {  
                 nxt[iz][ix]  = cur[iz][ix]*(aa[0][0]+aa[1][0])
                              + (cur[iz][ix-1]+cur[iz][ix+1])*aa[0][1]
                              + (cur[iz][ix-2]+cur[iz][ix+2])*aa[0][2]
                              + (cur[iz-1][ix]+cur[iz+1][ix])*aa[1][1]
                              + (cur[iz-2][ix]+cur[iz+2][ix])*aa[1][2];
             }
         }  
    
          
         nxt[isz+nb+2][isx+nb+2] += wav[it]/(v[isz+nb][isx+nb]*v[isz+nb][isx+nb]*dt*dt);

	 for (iz=2; iz < nzb-2; iz++) {  
             for (ix=2; ix < nxb-2; ix++) {
                 dercur[iz][ix]= derold[iz][ix] + nxt[iz][ix]*dt*v[iz-2][ix-2]*v[iz-2][ix-2];
                 nxt[iz][ix] = cur[iz][ix] + dercur[iz][ix]*dt; 
             }
         }

 
	 /*     nxt[isz+nb][isx+nb] += wav[it]; */
         
                 
	 for (iz=0; iz < nb; iz++) {  
             for (ix=nb+2; ix < nx+nb+2; ix++) {
                 nxt[iz+2][ix] *= wb[iz];
                 nxt[iz+nb+nz+2][ix] *= wb[nb-1-iz];
                 dercur[iz+2][ix] *= wb[iz];
                 dercur[iz+nb+nz+2][ix] *= wb[nb-1-iz];
             }
         }
	 for (iz=nb+2; iz < nz+nb+2; iz++) {  
             for (ix=0; ix < nb; ix++) {
                 nxt[iz][ix+2] *= wb[ix];
                 nxt[iz][ix+nx+nb+2] *= wb[nb-1-ix]; 
                 dercur[iz][ix+2] *= wb[ix];
                 dercur[iz][ix+nx+nb+2] *= wb[nb-1-ix];
             }
         }
	 for (iz=0; iz < nb; iz++) {  
             for (ix=0; ix < nb; ix++) {
                 nxt[iz+2][ix+2] *= wb[ix>iz?ix:iz];
                 dercur[iz+2][ix+2] *= wb[ix>iz?ix:iz];
                 nxt[iz+2][ix+nx+nb+2] *= wb[iz>(nb-1-ix)?iz:(nb-1-ix)]; 
                 dercur[iz+2][ix+nx+nb+2] *= wb[iz>(nb-1-ix)?iz:(nb-1-ix)]; 
                 nxt[iz+nb+nz+2][ix+2] *= wb[ix>(nb-1-iz)?ix:(nb-1-iz)];
                 dercur[iz+nb+nz+2][ix+2] *= wb[ix>(nb-1-iz)?ix:(nb-1-iz)];
                 nxt[iz+nb+nz+2][ix+nb+nz+2] *= wb[ix<iz?(nb-1-ix):(nb-1-iz)];
                 dercur[iz+nb+nz+2][ix+nb+nz+2] *= wb[ix<iz?(nb-1-ix):(nb-1-iz)];
             }
         }
	 for (iz=0; iz < nzb; iz++) {  
             for(ix=0; ix < nxb; ix++) {
	        old[iz][ix] = cur[iz][ix]; 
	        cur[iz][ix] = nxt[iz][ix]; 
	        derold[iz][ix] = dercur[iz][ix]; 
             }
         }
         for (iz=nb+2; iz<nz+nb+2; iz++){
             sf_floatwrite(nxt[iz]+nb+2,nx,out);
         }  
    }
    free(*aa);
    free(aa);
    free(*v);     
    free(*nxt);     
    free(*cur);     
    free(*old);     
    free(*dercur);     
    free(*derold);     
    free(v);     
    free(nxt);     
    free(cur);     
    free(old);     
    free(dercur);     
    free(derold);     
 
    exit(0); 
}           
           
