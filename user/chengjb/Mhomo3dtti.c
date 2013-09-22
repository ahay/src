/* Create a 3-D tilted TI model.

   Copyright (C) 2013 Tongji University, Shanghai, China 

   Authors: Jiubing Cheng
     
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

int main(int  argc,char **argv)
{
    int   ny,nx,nz;
    int   i,j,k;
    float dx,dy,dz,fx,fy,fz,vp,vs,ep,de,ga,th,ph;

    float ***v;

    sf_init(argc,argv);

    /* time samping paramter */
    if (!sf_getint("nx",&nx))   nx=101;
    if (!sf_getint("ny",&ny))   ny=101;
    if (!sf_getint("nz",&nz))   nz=101;
    if (!sf_getfloat("dx",&dx)) dx=0.008;
    if (!sf_getfloat("dy",&dy)) dy=0.008;
    if (!sf_getfloat("dz",&dz)) dz=0.008;
    if (!sf_getfloat("fx",&fx)) fx=0.0;
    if (!sf_getfloat("fy",&fy)) fy=0.0;
    if (!sf_getfloat("fz",&fz)) fz=0.0;
    if (!sf_getfloat("vp",&vp)) vp=3000.0;
    if (!sf_getfloat("vs",&vs)) vs=1200.0;
    if (!sf_getfloat("ep",&ep)) ep=0.2;
    if (!sf_getfloat("de",&de)) de=0.1;
    if (!sf_getfloat("ga",&ga)) ga=0.1;
    if (!sf_getfloat("th",&th)) th=0.0;
    if (!sf_getfloat("ph",&ph)) ph=0.0;

    sf_warning("nx= %d ny= %d nz= %d",nx,ny,nz);
    sf_warning("dx= %f dy= %f dz= %f",dx,dy,dz);
    sf_warning("vp0= %f vs0= %f eps= %f del=%f gam=%f the=%f pha=%f",vp,vs,ep,de,ga,th,ph);

    v = sf_floatalloc3(nz,nx,ny);

    /* setup I/O files */
    sf_file Fvp0, Fvs0, Fep, Fde, Fga;
    sf_file Fthe, Fphi;

    Fvp0 = sf_output ("out");  /* vp0 using standard output */
    Fvs0 = sf_output ("vs0");  /* vs0 */
    Fep = sf_output ("epsi");  /* epsi */
    Fde = sf_output ("del");  /* delta */
    Fga = sf_output ("gama");  /* gama */
    Fthe = sf_output ("the");  /* theta */
    Fphi = sf_output ("phi");  /* phai */

    puthead3x(Fvp0, nz, nx, ny, dz, dx, dy, fz, fx, fy);
    puthead3x(Fvs0, nz, nx, ny, dz, dx, dy, fz, fx, fy);
    puthead3x(Fep, nz, nx, ny, dz, dx, dy, fz, fx, fy);
    puthead3x(Fde, nz, nx, ny, dz, dx, dy, fz, fx, fy);
    puthead3x(Fga, nz, nx, ny, dz, dx, dy, fz, fx, fy);
    puthead3x(Fthe, nz, nx, ny, dz, dx, dy, fz, fx, fy);
    puthead3x(Fphi, nz, nx, ny, dz, dx, dy, fz, fx, fy);

    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          for(k=0;k<nz;k++) v[i][j][k] = vp;
          sf_floatwrite(v[i][j],nz,Fvp0);
        }

    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          for(k=0;k<nz;k++) v[i][j][k] = vs;
          sf_floatwrite(v[i][j],nz,Fvs0);
        }
    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          for(k=0;k<nz;k++) v[i][j][k] = ep;
          sf_floatwrite(v[i][j],nz,Fep);
        }
    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          for(k=0;k<nz;k++) v[i][j][k] = de;
          sf_floatwrite(v[i][j],nz,Fde);
        }
    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          for(k=0;k<nz;k++) v[i][j][k] = ga;
          sf_floatwrite(v[i][j],nz,Fga);
        }
    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          for(k=0;k<nz;k++) v[i][j][k] = th;
          sf_floatwrite(v[i][j],nz,Fthe);
        }
    for(i=0;i<ny;i++)
        for(j=0;j<nx;j++){
          for(k=0;k<nz;k++) v[i][j][k] = ph;
          sf_floatwrite(v[i][j],nz,Fphi);
        }
       
    sf_warning("create velocity model parameters ok");

    printf("ok3\n");

    free(**v);
		
    return 0;
}
