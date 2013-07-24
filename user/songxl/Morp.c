/* 2-D 10th-order Finite-difference dispersion*/
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
    int nx,  nz,  ny, k, aj;
    float dx, dz, dy, ox, oz, oy; /* dx2, dz2, dy2; */
    float vx, vy, vz, e1, e2, e3, wx, wy, wz, aa, bb, cc, r, mm;
    sf_file out, vel;
    float dkx, dkz, dky, dk2, dangle, gt, seta, phi;
    float kx, ky, kz, x, y, z;
    float pi=SF_PI, *ktmp, **vratio;
    float con2 = powf(2.0,1/3.0);

    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("in");   /* velocity */

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"o1",&oz)) oz=0.0;
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o2",&ox)) ox=0.0;
    if (!sf_histint(vel,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histfloat(vel,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(vel,"o3",&oy)) oy=0.0;
    if (!sf_getfloat("vx",&vx)) vx=2.1;
    if (!sf_getfloat("vy",&vy)) vy=2.05;
    if (!sf_getfloat("vz",&vz)) vz=2.0;
    if (!sf_getfloat("e1",&e1)) e1=0.3;
    if (!sf_getfloat("e2",&e2)) e2=0.1;
    if (!sf_getfloat("e3",&e3)) e3=1.0;
    if (!sf_getfloat("phi",&phi)) phi=45.0;
        phi *= pi/180.0;

    dkx = 1.0/(dx*nx);
    dkz = 1.0/(dz*nz);
    dky = 1.0/(dy*ny);
    dk2=sqrtf(dkx*dkx+dkz*dkz+dky*dky)*2.0*pi;
    dangle = pi/4;
/*    phi = pi/4; */
    sf_putint(out,"n1",nz);
    sf_putfloat(out,"d1",dk2);
    sf_putfloat(out,"o1",oz); 
    sf_putint(out,"n2",3);
    sf_putfloat(out,"d2",dangle);
    sf_putint(out,"n3",1);

    vratio =  sf_floatalloc2(nz,3);
    ktmp =  sf_floatalloc(nz);
    for (k=0; k < nz; k++) {
        ktmp[k] = k*dk2;
    }


/*    dx2 = dx*dx;
      dz2 = dz*dz;
      dy2 = dy*dy; */
    wx = vx*vx;
    wy = vy*vy;
    wz = vz*vz;
    for (aj=0; aj<3; aj++) {
        seta = dangle*aj;
        for (k=1; k<nz; k++) {
            gt = 0;
            kz = ktmp[k]*cosf(seta);
            kx = ktmp[k]*sinf(seta)*cosf(phi);
            ky = ktmp[k]*sinf(seta)*sinf(phi);
             x = kx*kx;
             y = ky*ky;
             z = kz*kz;
            aa = (2*e1+1)*wx*x+(2*e2+1)*wy*y+wz*z;
            bb = wx*wx*x*y*(2*e1*e3+e3)*(2*e1*e3+e3)-wx*wy*(2*e1+1)*(2*e2+1)*x*y-2*wx*wz*e1*x*z-2*wy*wz*e2*y*z;
            cc = (wz*z)*(wx*x)*y*(-(wx)*(2*e1*e3+e3)*(2*e1*e3+e3)+2*(vx*vy)*e3*(2*e1+1)-(wy)*(1-4*e1*e2));
             r = (81*cc+6*aa*(2*aa*aa+9*bb))*cc-3*bb*bb*(aa*aa+4*bb);
             r = sqrtf(fabsf(r))-9*cc;
            mm = -2*aa*aa*aa+3*r-9*aa*bb;
            if (mm<0) r = -powf(-mm,(1.0/3.0));
               else r = powf(mm,(1.0/3.0));
            if (fabsf(r) < 0.000001) {r = 0.0;}
               else { r = 1/6.0*(-con2*con2*r-2*con2*(aa*aa+3*bb)/r+2*aa);}
            r = sqrtf(fabsf(r));
            gt = r/ktmp[k];
            vratio[aj][k]=gt;
        }
        vratio[aj][0]=vratio[aj][1];
    }
    sf_floatwrite(vratio[0],nz*3,out);
    exit(0); 
}           
           
