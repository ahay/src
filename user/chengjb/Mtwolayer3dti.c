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
    int    nx, ny, nz;
    float  dx, dy, dz; 
    float  vp0_1, vs0_1, eps_1, del_1, gam_1, the_1, phi_1;
    float  vp0_2, vs0_2, eps_2, del_2, gam_2, the_2, phi_2;
    float  ***w;
    int   i,j,k;
    int nnz;
    /* setup I/O files */
    sf_file Fo1, Fo2, Fo3, Fo4, Fo5, Fo6, Fo7;

    sf_init(argc,argv);

    if (!sf_getint("ny",&ny)) ny=201;
    if (!sf_getint("nx",&nx)) nx=201;
    if (!sf_getint("nz",&nz)) nz=201;
    if (!sf_getfloat("dy",&dy)) dy=0.008;
    if (!sf_getfloat("dx",&dx)) dx=0.008;
    if (!sf_getfloat("dz",&dz)) dz=0.008;
    if (!sf_getfloat("vp0_1",&vp0_1)) vp0_1=3000.0;
    if (!sf_getfloat("vs0_1",&vs0_1)) vs0_1=1200.0;
    if (!sf_getfloat("eps_1",&eps_1)) eps_1=0.2;
    if (!sf_getfloat("del_1",&del_1)) del_1=0.1;
    if (!sf_getfloat("gam_1",&gam_1)) gam_1=0.0;
    if (!sf_getfloat("the_1",&the_1)) the_1=0.0;
    if (!sf_getfloat("phi_1",&phi_1)) phi_1=0.0;
    if (!sf_getfloat("vp0_2",&vp0_2)) vp0_2=3000.0;
    if (!sf_getfloat("vs0_2",&vs0_2)) vs0_2=1200.0;
    if (!sf_getfloat("eps_2",&eps_2)) eps_2=0.2;
    if (!sf_getfloat("del_2",&del_2)) del_2=0.1;
    if (!sf_getfloat("gam_2",&gam_2)) gam_2=0.0;
    if (!sf_getfloat("the_2",&the_2)) the_2=30.0;  /* Unit: degree */
    if (!sf_getfloat("phi_2",&phi_2)) phi_2=0.0;

    sf_warning("nx= %d ny=%d nz= %d",nx,ny,nz);
    sf_warning("dx= %f dy=%f dz= %f",dx,dy,dz);
    sf_warning("vp0_1= %f vs0_1= %f eps_1= %f del_1=%f gam_1=%f the_1=%f phi_1=%f",vp0_1,vs0_1,eps_1,del_1,gam_1,the_1,phi_1);
    sf_warning("vp0_2= %f vs0_2= %f eps_2= %f del_2=%f gam_2=%f the_2=%f phi_1=%f",vp0_2,vs0_2,eps_2,del_2,gam_2,the_2,phi_2);

    w = sf_floatalloc3(nz,nx, ny);

    zero3float(w,nz,nx,ny);

    Fo1 = sf_output("out"); 
    Fo2 = sf_output("vs0"); 
    Fo3 = sf_output("epsi");
    Fo4 = sf_output("del"); 
    Fo5 = sf_output("gam"); 
    Fo6 = sf_output("the"); 
    Fo7 = sf_output("phi"); 

    puthead3x(Fo1, nz, nx, ny, dz, dx, dy, 0.0, 0.0, 0.0);
    puthead3x(Fo2, nz, nx, ny, dz, dx, dy, 0.0, 0.0, 0.0);
    puthead3x(Fo3, nz, nx, ny, dz, dx, dy, 0.0, 0.0, 0.0);
    puthead3x(Fo4, nz, nx, ny, dz, dx, dy, 0.0, 0.0, 0.0);
    puthead3x(Fo5, nz, nx, ny, dz, dx, dy, 0.0, 0.0, 0.0);
    puthead3x(Fo6, nz, nx, ny, dz, dx, dy, 0.0, 0.0, 0.0);
    puthead3x(Fo7, nz, nx, ny, dz, dx, dy, 0.0, 0.0, 0.0);

    nnz=nz*7/12;

    /* Vp0 */
    for(k=0;k<ny;k++)
        for(i=0;i<nx;i++)
        {
	    for(j=0;j<nnz;j++)
		w[k][i][j] = vp0_1;
	    for(j=nnz;j<nz;j++)
		w[k][i][j] = vp0_2;
           
	    sf_floatwrite(w[k][i], nz, Fo1);
        }
    /* Vs0 */
    for(k=0;k<ny;k++)
        for(i=0;i<nx;i++)
        {
	    for(j=0;j<nnz;j++)
		w[k][i][j] = vs0_1;
	    for(j=nnz;j<nz;j++)
		w[k][i][j] = vs0_2;
           
	    sf_floatwrite(w[k][i], nz, Fo2);
        }
    /* Epsi */
    for(k=0;k<ny;k++)
        for(i=0;i<nx;i++)
        {
	    for(j=0;j<nnz;j++)
		w[k][i][j] = eps_1;
	    for(j=nnz;j<nz;j++)
		w[k][i][j] = eps_2;
           
	    sf_floatwrite(w[k][i], nz, Fo3);
        }
    /* Delta */ 
    for(k=0;k<ny;k++)
        for(i=0;i<nx;i++)
        {
	    for(j=0;j<nnz;j++)
		w[k][i][j] = del_1;
	    for(j=nnz;j<nz;j++)
		w[k][i][j] = del_2;
           
	    sf_floatwrite(w[k][i], nz, Fo4);
        }
    /*  Gama */
    for(k=0;k<ny;k++)
        for(i=0;i<nx;i++)
        {
	    for(j=0;j<nnz;j++)
		w[k][i][j] = gam_1;
	    for(j=nnz;j<nz;j++)
		w[k][i][j] = gam_2;
           
	    sf_floatwrite(w[k][i], nz, Fo5);
        }
    /*  Theta */
    for(k=0;k<ny;k++)
        for(i=0;i<nx;i++)
        {
	    for(j=0;j<nnz;j++)
		w[k][i][j] = the_1;
	    for(j=nnz;j<nz;j++)
		w[k][i][j] = the_2;
           
	    sf_floatwrite(w[k][i], nz, Fo6);
        }
    /* Phai */ 
    for(k=0;k<ny;k++)
        for(i=0;i<nx;i++)
        {
	    for(j=0;j<nnz;j++)
		w[k][i][j] = phi_1;
	    for(j=nnz;j<nz;j++)
		w[k][i][j] = phi_2;
           
	    sf_floatwrite(w[k][i], nz, Fo7);
        }
    free(**w);

    exit(0);
}
