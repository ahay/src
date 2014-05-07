/* 2-D two-components wavefield modeling using pseudo-pure mode P-wave equation in VTI media.
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng and Wei Kang
     
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

/* calculate phase, group velocity and angle  */
#include "zero.h"
#include "puthead.h"

int main(int argc, char* argv[])
{

        int    nx, nz;
        float  dx, dz; 
        float  vp0_1, vs0_1, eps_1, del_1, the_1;
        float  vp0_2, vs0_2, eps_2, del_2, the_2;
        float  vp0_3, vs0_3, eps_3, del_3, the_3;
        float  **w;
	int   i,j;
        sf_file Fo1, Fo2, Fo3, Fo4, Fo5, Fo6;
	int nx1, nx2;

        sf_init(argc,argv);

        if (!sf_getint("nx",&nx)) nx=301;
        if (!sf_getint("nz",&nz)) nz=301;
        if (!sf_getfloat("dx",&dx)) dx=0.008;
        if (!sf_getfloat("dz",&dz)) dz=0.008;
        if (!sf_getfloat("vp0_1",&vp0_1)) vp0_1=3000.0;
        if (!sf_getfloat("vs0_1",&vs0_1)) vs0_1=1200.0;
        if (!sf_getfloat("eps_1",&eps_1)) eps_1=0.2;
        if (!sf_getfloat("del_1",&del_1)) del_1=0.1;
        if (!sf_getfloat("the_1",&the_1)) the_1=0.0;
        if (!sf_getfloat("vp0_2",&vp0_2)) vp0_2=3000.0;
        if (!sf_getfloat("vs0_2",&vs0_2)) vs0_2=1200.0;
        if (!sf_getfloat("eps_2",&eps_2)) eps_2=0.2;
        if (!sf_getfloat("del_2",&del_2)) del_2=0.1;
        if (!sf_getfloat("the_2",&the_2)) the_2=30.0;  /* Unit: degree */
        if (!sf_getfloat("vp0_3",&vp0_3)) vp0_3=3000.0;
        if (!sf_getfloat("vs0_3",&vs0_3)) vs0_3=1200.0;
        if (!sf_getfloat("eps_3",&eps_3)) eps_3=0.2;
        if (!sf_getfloat("del_3",&del_3)) del_3=0.1;
        if (!sf_getfloat("the_3",&the_3)) the_3=0.0;  /* Unit: degree */

        sf_warning("nx= %d nz= %d",nx,nz);
        sf_warning("dx= %f dz= %f",dx,dz);
        sf_warning("vp0_1= %f vs0_1= %f eps_1= %f del_1=%f the_1=%f",vp0_1,vs0_1,eps_1,del_1,the_1);
        sf_warning("vp0_2= %f vs0_2= %f eps_2= %f del_2=%f the_2=%f",vp0_2,vs0_2,eps_2,del_2,the_2);
        sf_warning("vp0_3= %f vs0_3= %f eps_3= %f del_3=%f the_3=%f",vp0_3,vs0_3,eps_3,del_3,the_3);

        w = sf_floatalloc2(nz,nx);

	zero2float(w,nz,nx);

        /* setup I/O files */
        Fo1 = sf_output("out"); /* Vp0 */
        Fo2 = sf_output("vs0"); /* Vs0 */
        Fo3 = sf_output("epsi"); /* Epsilon */
        Fo4 = sf_output("del"); /* Delta */
        Fo5 = sf_output("the"); /* theta */
        Fo6 = sf_output("interf"); /* interface */

        puthead2(Fo1, nz, nx, dz, 0.0, dx, 0.0);
        puthead2(Fo2, nz, nx, dz, 0.0, dx, 0.0);
        puthead2(Fo3, nz, nx, dz, 0.0, dx, 0.0);
        puthead2(Fo4, nz, nx, dz, 0.0, dx, 0.0);
        puthead2(Fo5, nz, nx, dz, 0.0, dx, 0.0);
        puthead2(Fo6, nz, nx, dz, 0.0, dx, 0.0);

        nx1=nx/3;
        nx2=(nx*2)/3;

        /* Vp0 */
        for(i=0;i<nz;i++)
        {
           for(j=0;j<nx1;j++)
              w[j][i] = vp0_1;
           for(j=nx1;j<nx2;j++)
              w[j][i] = vp0_2;
           for(j=nx2;j<nx;j++)
              w[j][i] = vp0_3;
        }
        for(j=0;j<nx;j++)
          sf_floatwrite(w[j], nz, Fo1);

        /* Vs0 */
        for(i=0;i<nz;i++)
        {
           for(j=0;j<nx1;j++)
              w[j][i] = vs0_1;
           for(j=nx1;j<nx2;j++)
              w[j][i] = vs0_2;
           for(j=nx2;j<nx;j++)
              w[j][i] = vs0_3;
        }
        for(j=0;j<nx;j++)
           sf_floatwrite(w[j], nz, Fo2);

        /* Epsi */
        for(i=0;i<nz;i++)
        {
           for(j=0;j<nx1;j++)
              w[j][i] = eps_1;
           for(j=nx1;j<nx2;j++)
              w[j][i] = eps_2;
           for(j=nx2;j<nx;j++)
              w[j][i] = eps_3;
        }
        for(j=0;j<nx;j++)
           sf_floatwrite(w[j], nz, Fo3);

        /* Delta */ 
        for(i=0;i<nz;i++)
        {
           for(j=0;j<nx1;j++)
              w[j][i] = del_1;
           for(j=nx1;j<nx2;j++)
              w[j][i] = del_2;
           for(j=nx2;j<nx;j++)
              w[j][i] = del_2;
        }
        for(j=0;j<nx;j++)
           sf_floatwrite(w[j], nz, Fo4);

        /*  Theta */
        for(i=0;i<nz;i++)
        {
           for(j=0;j<nx1;j++)
              w[j][i] = the_1;
           for(j=nx1;j<nx2;j++)
              w[j][i] = the_2;
           for(j=nx2;j<nx;j++)
              w[j][i] = the_2;
        }
        for(j=0;j<nx;j++)
           sf_floatwrite(w[j], nz, Fo5);

        /*  Interface */
        for(i=0;i<nz;i++)
        {
           for(j=0;j<nx;j++)
              w[j][i] = 0.0;

           w[nx1][i] = 1;
           w[nx2][i] = 1;
        }
        for(j=0;j<nx;j++)
           sf_floatwrite(w[j], nz, Fo6);

        free(*w);

	exit(0);
}
