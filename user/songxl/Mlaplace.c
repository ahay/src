/* 2-D finite-difference Laplacian */
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
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx,  nz, ix, iz, nbt, nbb, nbl, nbr, nxb, nzb;
    float dx, dz;
    float **new,  **old;
    sf_file out, img;
     

    sf_init(argc,argv);
    out = sf_output("out");
    img= sf_input("in");   /* input file*/

    if (SF_FLOAT != sf_gettype(img)) sf_error("Need float input");
    if (!sf_histint(img,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(img,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(img,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(img,"d2",&dz)) sf_error("No d2= in input");

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);

    nbl =2;
    nbr =2;
    nbt =2;
    nbb =2; 
    nxb = nx + nbl + nbr;
    nzb = nz + nbt + nbb;




    old    =  sf_floatalloc2(nxb,nzb);
    new    =  sf_floatalloc2(nx,nz);
    
    bd_init(nx,nz,nbt,nbb,nbl,nbr,0,0,0,0);

    /*input & extend image*/
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(old[iz]+nbl,nx,img);
         for (ix=0; ix<nbl; ix++){
             old[iz][ix] = old[iz][nbl];
         }
         for (ix=0; ix<nbr; ix++){
             old[iz][nx+nbl+ix] = old[iz][nx+nbl-1];
         }     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            old[iz][ix] = old[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            old[nz+nbt+iz][ix] = old[nz+nbt-1][ix];
        }
    }

 
    for (iz=0; iz < nz; iz++){
         for (ix=0; ix < nx; ix++) {
             new[iz][ix] = (old[iz+nbt][ix+nbl-2]-16.0*old[iz+nbt][ix+nbl-1]+30.0*old[iz+nbt][ix+nbl]-16.0*old[iz+nbt][ix+nbl+1]+old[iz+nbt][ix+nbl+2])/(12.0*dx*dx)
                        + (old[iz+nbt-2][ix+nbl]-16.0*old[iz+nbt-1][ix+nbl]+30.0*old[iz+nbt][ix+nbl]-16.0*old[iz+nbt+1][ix+nbl]+old[iz+nbt+2][ix+nbl])/(12.0*dz*dz);
         }
    }

    sf_floatwrite(new[0],nx*nz,out);
    bd_close();
    free(*new);     
    free(*old);     
    free(new);     
    free(old);     
    exit(0); 
}           
           
