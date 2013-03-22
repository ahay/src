/* Correct projection deviation in X-domian for 3-D pseudo-pure P-wave field in homogeneous ORT media.

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
#include "zero.h"

/* spectra multiply for wave-mode separation*/
#include "kykxkztaper.h"
#include "spectramultipy.h"

int main(int  argc,char **argv)
{
    sf_init(argc,argv);

    sf_file Fix, Fiy, Fiz;
    sf_file Fi1, Fi2, Fi3;
    sf_file Fo1, Fo2;

    clock_t t1, t2;
    float   timespent;

    t1=clock();

   int isep;

   if (!sf_getint("isep",&isep)) isep=1;

    /* setup I/O files */
    Fix = sf_input("in");
    Fiy = sf_input("apvyy");
    Fiz = sf_input("apvzz");

    Fi1 = sf_input("PseudoPurePx"); /* pseudo-pure P-wave x-component */
    Fi2 = sf_input("PseudoPurePy"); /* pseudo-pure P-wave y-component */
    Fi3 = sf_input("PseudoPurePz"); /* pseudo-pure P-wave z-component */
    Fo1 = sf_output("out");               /* pseudo-pure scalar P-wave */
    Fo2 = sf_output("PseudoPureSepP");    /* separated scalar P-wave */

    int     nx, ny, nz, nxf, nyf, nzf, hnx, hny, hnz, i, j, k;
    float   fx, fy, fz, dx, dy, dz, dxf, dyf, dzf;

    /* Read/Write axes */
    sf_axis az, ax, ay;

    az = sf_iaxa(Fi1,1); nz = sf_n(az); dz = sf_d(az)*1000.0;
    ax = sf_iaxa(Fi1,2); nx = sf_n(ax); dx = sf_d(ax)*1000.0;
    ay = sf_iaxa(Fi1,3); ny = sf_n(ay); dy = sf_d(ay)*1000.0;
    fy=sf_o(ay)*1000.0;
    fx=sf_o(ax)*1000.0;
    fz=sf_o(az)*1000.0;

    sf_warning("nx=%d ny=%d nz=%d dx=%f dy=%f dz=%f",nx,ny,nz,dx,dy,dz);

    puthead3x(Fo1, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);
    puthead3x(Fo2, nz, nx, ny, dz/1000.0, dx/1000.0, dy/1000.0, 0.0, 0.0, 0.0);

    float*** px=sf_floatalloc3(nz,nx,ny);
    float*** py=sf_floatalloc3(nz,nx,ny);
    float*** pz=sf_floatalloc3(nz,nx,ny);
    float*** p=sf_floatalloc3(nz,nx,ny);

    int iy, ix, iz;

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

if(isep==1){
    /* Read axes */
    sf_axis azf, axf, ayf;
    azf = sf_iaxa(Fix,1); nzf = sf_n(azf); dzf = sf_d(azf)*1000.0;
    axf = sf_iaxa(Fix,2); nxf = sf_n(axf); dxf = sf_d(axf)*1000.0;
    ayf = sf_iaxa(Fix,3); nyf = sf_n(ayf); dyf = sf_d(ayf)*1000.0;

    if(dx!=dxf){
       sf_warning("dx= %f dxf=%f",dx,dxf);
       sf_warning("filter and data spatial sampling don't match in x-axis");
       exit(0);
    }
    if(dy!=dyf){
       sf_warning("dy= %f dyf=%f",dy,dyf);
       sf_warning("filter and data spatial sampling don't match in y-axis");
       exit(0);
    }
    if(dz!=dzf){
       sf_warning("dz= %f dzf=%f",dz,dzf);
       sf_warning("filter and data spatial sampling don't match in z-axis");
       exit(0);
    }

    float*** apvxx=sf_floatalloc3(nzf,nxf,nyf);
    float*** apvyy=sf_floatalloc3(nzf,nxf,nyf);
    float*** apvzz=sf_floatalloc3(nzf,nxf,nyf);

    hnx=nxf/2;
    hny=nyf/2;
    hnz=nzf/2;

    sf_warning("nxf=%d nyf=%d nzf=%d dxf=%f dyf=%f dzf=%f",nxf,nyf,nzf,dxf,dyf,dzf);

    for(iy=0;iy<nyf;iy++)
    for(ix=0;ix<nxf;ix++)
    {
      sf_floatread(apvxx[iy][ix],nzf,Fix);
      sf_floatread(apvyy[iy][ix],nzf,Fiy);
      sf_floatread(apvzz[iy][ix],nzf,Fiz);
    }

    float*** pxc=sf_floatalloc3(nz,nx,ny);
    float*** pyc=sf_floatalloc3(nz,nx,ny);
    float*** pzc=sf_floatalloc3(nz,nx,ny);
    zero3float(pxc, nz,nx,ny);
    zero3float(pyc, nz,nx,ny);
    zero3float(pzc, nz,nx,ny);

    int l, g, h, ll, gg, hh;

    for(j=0;j<ny;j++){
      sf_warning("ny=%d iy=%d",ny,j);
      for(i=0;i<nx;i++)
      for(k=0;k<nz;k++)
      {

        for(g=-hny; g<=hny; g++){
            gg=g+hny;
            for(l=-hnx; l<=hnx; l++){
                ll=l+hnx;
                for(h=-hnz; h<=hnz; h++)
                {
                    hh=h+hnz;
                    if(i+l>=0 && i+l<nx && j+g>=0 && j+g<ny && k+h>=0 && k+h<nz){
                        pyc[i][j][k]+=py[i+l][j+g][k+h]*apvxx[gg][ll][hh];
                        pxc[i][j][k]+=px[i+l][j+g][k+h]*apvyy[gg][ll][hh];
                        pzc[i][j][k]+=pz[i+l][j+g][k+h]*apvzz[gg][ll][hh];
                    }
                } // h loop
             } // g loop
          }// l oop
       }
    }
/*
    int iix, iiy, iiz;

    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
       for(iz=0;iz<nz;iz++){
          pxc[iy][ix][iz]=px[iy][ix][iz];
          pyc[iy][ix][iz]=py[iy][ix][iz];
          pzc[iy][ix][iz]=pz[iy][ix][iz];
          iix=ix-nx/2;
          iiy=iy-ny/2;
          iiz=iz-nz/2;
          if(sqrt(1.0*(iix*iix+iiy*iiy+iiz*iiz))<3.0*nx/8)
          {
             pxc[iy][ix][iz] = 0.0;
             pyc[iy][ix][iz] = 0.0;
             pzc[iy][ix][iz] = 0.0;
          }
       }
*/
    for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++){
       for(iz=0;iz<nz;iz++)
          p[iy][ix][iz] = pxc[iy][ix][iz] + pyc[iy][ix][iz] + pzc[iy][ix][iz];
       sf_floatwrite(p[iy][ix],nz,Fo2);
    }
    free(**pxc);
    free(**pyc);
    free(**pzc);
    free(**apvxx);
    free(**apvyy);
    free(**apvzz);
}// endif
    t2=clock();
    timespent=(float)(t2-t1)/CLOCKS_PER_SEC;
    sf_warning("Computation time (Filtering): %f (second)",timespent);

    free(**px);
    free(**py);
    free(**pz);

    free(**p);
		
    return 0;
}
