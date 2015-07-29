/* Correct projection deviation in K-domian for 3-D pseudo-pure P-wave field in homogeneous ORT media.

   Refernces:
             Cheng et al. (15th IWSA, 2012);
             Cheng and Kang (SEG Abstract, 2012);
             Kang and Cheng (SEG Abstract, 2012)
             Wang et al.(SEG Abstract, 2012)      

   Copyright (C) 2012 Tongji University, Shanghai, China 

   Authors: Jiubing Cheng, Tengfei Wang and Wei Kang
     
   This code is first written by Tengfei Wang at Tongji University,
   and then optimzied by Jiubing Cheng for Madagascar version at BEG,
   University of Texas at Austin.

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

/* prepared head files by myself */
#include "_cjb.h"

/* head files aumatically produced from *.c */
#include "puthead.h"

/* spectra multiply for wave-mode separation*/
#include "kykxkztaper.h"
#include "spectramultipy.h"

int main(int  argc,char **argv)
{
    sf_file Fix, Fiy, Fiz;
    sf_file Fi1, Fi2, Fi3;
    sf_file Fo1, Fo2;

    clock_t t1, t2;
    float   timespent;
    int     nx, ny, nz, nkx, nky, nkz;
    float   dx, dy, dz;

    /* Read/Write axes */
    sf_axis az, ax, ay;

    float*** px;
    float*** py;
    float*** pz;
    float*** p;

    int iy, ix, iz, iflag=1;
    sf_axis akz, akx, aky;
    int *ijkx, *ijky, *ijkz;

    float*** apvx;
    float*** apvy;
    float*** apvz;

    t1=clock();

    /* time samping paramter */
    /* setup I/O files */
    sf_init(argc,argv);

    Fix = sf_input("in");
    Fiy = sf_input("apvy");
    Fiz = sf_input("apvz");

    Fi1 = sf_input("PseudoPurePx"); /* pseudo-pure P-wave x-component */
    Fi2 = sf_input("PseudoPurePy"); /* pseudo-pure P-wave y-component */
    Fi3 = sf_input("PseudoPurePz"); /* pseudo-pure P-wave z-component */
    Fo1 = sf_output("out");               /* pseudo-pure scalar P-wave */
    Fo2 = sf_output("PseudoPureSepP");    /* separated scalar P-wave */



    az = sf_iaxa(Fi1,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
    ax = sf_iaxa(Fi1,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
    ay = sf_iaxa(Fi1,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
/*
    fy=sf_o(ay)*1000.0;
    fx=sf_o(ax)*1000.0;
    fz=sf_o(az)*1000.0;
*/

    sf_warning("nx=%d ny=%d nz=%d dx=%f dy=%f dz=%f",nx,ny,nz,dx,dy,dz);

    puthead3x(Fo1, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
    puthead3x(Fo2, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);

    px=sf_floatalloc3(nz,nx,ny);
    py=sf_floatalloc3(nz,nx,ny);
    pz=sf_floatalloc3(nz,nx,ny);
    p=sf_floatalloc3(nz,nx,ny);

    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
    {
      sf_floatread(px[iy][ix],nz,Fi1);
      sf_floatread(py[iy][ix],nz,Fi2);
      sf_floatread(pz[iy][ix],nz,Fi3);
    }

    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++){
       for(iz=0;iz<nz;iz++)
          p[iy][ix][iz] = px[iy][ix][iz] + py[iy][ix][iz] + pz[iy][ix][iz];
       sf_floatwrite(p[iy][ix],nz,Fo1);
    }

    /* Read axes */
    akz = sf_iaxa(Fix,1); nkz = sf_n(akz); /* dkz = sf_d(akz); */
    akx = sf_iaxa(Fix,2); nkx = sf_n(akx); /* dkx = sf_d(akx); */
    aky = sf_iaxa(Fix,3); nky = sf_n(aky); /* dky = sf_d(aky); */

    if(nx!=nkx){
       sf_warning("nx= %d nkx=%d",nx,nkx);
       sf_warning("filter and data length don't match in x-axis");
       exit(0);
    }
    if(ny!=nky){
       sf_warning("ny= %d nky=%d",ny,nky);
       sf_warning("filter and data length don't match in y-axis");
       exit(0);
    }
    if(nz!=nkz){
       sf_warning("nz= %d nkz=%d",nz,nkz);
       sf_warning("filter and data length don't match in z-axis");
       exit(0);
    }

    ijkx = sf_intalloc(nx);
    ijky = sf_intalloc(ny);
    ijkz = sf_intalloc(nz);

    ikxikyikz(ijkx, ijky, ijkz, nx, ny, nz);

    apvx=sf_floatalloc3(nz,nx,ny);
    apvy=sf_floatalloc3(nz,nx,ny);
    apvz=sf_floatalloc3(nz,nx,ny);

    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
    {
      sf_floatread(apvx[iy][ix],nz,Fix);
      sf_floatread(apvy[iy][ix],nz,Fiy);
      sf_floatread(apvz[iy][ix],nz,Fiz);
    }
    spec3dmultiply(px, apvx, nx, ny, nz, ijkx, ijky, ijkz, iflag);
    spec3dmultiply(py, apvy, nx, ny, nz, ijkx, ijky, ijkz, iflag);
    spec3dmultiply(pz, apvz, nx, ny, nz, ijkx, ijky, ijkz, iflag);
/*
    int iix, iiy, iiz;

    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
       for(iz=0;iz<nz;iz++){
          iix=ix-nx/2;
          iiy=iy-ny/2;
          iiz=iz-nz/2;
          if(sqrt(1.0*(iix*iix+iiy*iiy+iiz*iiz))<3.0*nx/8)
          {
             px[iy][ix][iz] = 0.0;            
             py[iy][ix][iz] = 0.0;            
             pz[iy][ix][iz] = 0.0;            
          }
       }
*/
      
    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++){
       for(iz=0;iz<nz;iz++)
          p[iy][ix][iz] = px[iy][ix][iz] + py[iy][ix][iz] + pz[iy][ix][iz];
       sf_floatwrite(p[iy][ix],nz,Fo2);
    }

    t2=clock();
    timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
    sf_warning("Computation time (Filtering): %f (second)",timespent);

    free(**apvx);
    free(**apvy);
    free(**apvz);
    free(**px);
    free(**py);
    free(**pz);
    free(**p);
		
    free(ijkx);
    free(ijky);
    free(ijkz);

    exit(0);
}
