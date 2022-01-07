/* Modeling of pure acoustic wave in 2-D transversely isotropic meida using optimized pseudo-Laplacian operator */
/*
  Copyright (C)2014 Institute of Geology and Geophysics, Chinese Academy of Sciences (Jun Yan) 
					2009 University of Texas at Austin
  
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
#include "abcpass2.h"
#include "opostep2.h"
#include "fft2.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nz, nt, ix, iz, it, nbt, nbb, nxl, nxr,  nxb, nzb, isx, isz;
    float dt, dx, dz, o1, o2;
    float **old,  **cur,  **tmp, *wav;
    float  **v, **vtmp, v0, **sigma, **delta, **seta;
    float ***aa, w, g1, g2, czt, czb, cxl, cxr; /* top, bottom, left, right */
    float ax, az, factor;
    sf_file out, vel, source, fsigma, fdelta, fseta, Gxxl, Gxxr, Gzzl, Gzzr,  Gxxxxl, Gxxxxr, Gzzzzl, Gzzzzr, Gxzzzl, Gxzzzr, Gxxxzl, Gxxxzr, Gxxzzl, Gxxzzr;
    float **Cxxl, **Cxxr, **Czzl, **Czzr,  **Cxxxxl, **Cxxxxr, **Czzzzl, **Czzzzr, **Cxzzzl, **Cxzzzr, **Cxxxzl, **Cxxxzr, **Cxxzzl, **Cxxzzr;
    int opt, snap, nsnap;    /* optimal padding */
	int nkxz,nkxx,nkzz,nxzb,nxzb2,n2,m2;
	int pad1;
	bool cmplx;
     
    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    fsigma = sf_input("sigma");   /* velocity */
    fdelta = sf_input("delta");   /* velocity */
    fseta  = sf_input("seta");   /* velocity */
    source = sf_input("in");   /* source wavlet*/
    
	Gxxl  = sf_input("Gxxl");
	Gxxr  = sf_input("Gxxr");
	Gzzl  = sf_input("Gzzl");
	Gzzr  = sf_input("Gzzr");
	Gxxxxl  = sf_input("Gxxxxl");
	Gxxxxr  = sf_input("Gxxxxr");
	Gzzzzl  = sf_input("Gzzzzl");
	Gzzzzr  = sf_input("Gzzzzr");
	Gxzzzl  = sf_input("Gxzzzl");
	Gxzzzr  = sf_input("Gxzzzr");
	Gxxxzl  = sf_input("Gxxxzl");
	Gxxxzr  = sf_input("Gxxxzr");
	Gxxzzl  = sf_input("Gxxzzl");
	Gxxzzr  = sf_input("Gxxzzr");

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");

	if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(vel,"o2",&o2)) o2=0.0;
    /*  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input"); */
    /*  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input"); */
    
	if (!sf_getint("opt",&opt)) opt=1;
    /* if y, determine optimal size for efficiency */

    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    if (!sf_getint("nbt",&nbt)) nbt=30;
    if (!sf_getint("nbb",&nbb)) nbb=30;
    if (!sf_getint("nxl",&nxl)) nxl=30;
    if (!sf_getint("nxr",&nxr)) nxr=30;
	
	/* assume ABC pars are the same */
	if (nbt != nbb || nxl != nxr || nbt!=nxl) 
		sf_error("ABC pars are not the same");

    if (!sf_getfloat("czt",&czt))  czt = 0.01; /*decaying parameter*/
    if (!sf_getfloat("czb",&czb))  czb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxl",&cxl)) cxl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxr",&cxr)) cxr = 0.01; /*decaying parameter*/

    if (!sf_getfloat("ax",&ax)) ax= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("az",&az)) az= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("factor",&factor)) factor= 5.0/6.0; /*suppress HF parameter*/

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

	 if (!sf_getint("snap",&snap)) snap=1;
	 nsnap=0;
	 for (it=0; it < nt; it++) {
			if (it%snap == 0) nsnap++;
	  }

    sf_putfloat(out,"d1",dz);
    sf_putfloat(out,"d2",dx);
    sf_putfloat(out,"d3",dt*snap);
    sf_putfloat(out,"o1",o1); 
    sf_putfloat(out,"o2",o2); 
    sf_putfloat(out,"o3",0.0); 

    nxb = nx + nxl + nxr;
    nzb = nz + nbt + nbb;
    
	 sf_putint(out,"n1",nzb);
    sf_putint(out,"n2",nxb);
    sf_putint(out,"n3",nsnap);

	 nkxz=fft2_init(cmplx,pad1,nzb,nxb,&nkzz,&nkxx);
	 nxzb = nxb*nzb;
	 nxzb2 = nkxx * nkzz;


    if (!sf_histint(Gxxl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Gxxl,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(Gxxr,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(Gxxr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
 
    Cxxl = sf_floatalloc2(nxzb,m2);
    Cxxr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Cxxl[0],nxzb*m2,Gxxl);
    sf_floatread(Cxxr[0],m2*nkxz,Gxxr);

    Czzl = sf_floatalloc2(nxzb,m2);
    Czzr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Czzl[0],nxzb*m2,Gzzl);
    sf_floatread(Czzr[0],m2*nkxz,Gzzr);


    Cxxxxl = sf_floatalloc2(nxzb,m2);
    Cxxxxr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Cxxxxl[0],nxzb*m2,Gxxxxl);
    sf_floatread(Cxxxxr[0],m2*nkxz,Gxxxxr);

    Czzzzl = sf_floatalloc2(nxzb,m2);
    Czzzzr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Czzzzl[0],nxzb*m2,Gzzzzl);
    sf_floatread(Czzzzr[0],m2*nkxz,Gzzzzr);


    Cxzzzl = sf_floatalloc2(nxzb,m2);
    Cxzzzr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Cxzzzl[0],nxzb*m2,Gxzzzl);
    sf_floatread(Cxzzzr[0],m2*nkxz,Gxzzzr);

    Cxxxzl = sf_floatalloc2(nxzb,m2);
    Cxxxzr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Cxxxzl[0],nxzb*m2,Gxxxzl);
    sf_floatread(Cxxxzr[0],m2*nkxz,Gxxxzr);


    Cxxzzl = sf_floatalloc2(nxzb,m2);
    Cxxzzr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Cxxzzl[0],nxzb*m2,Gxxzzl);
    sf_floatread(Cxxzzr[0],m2*nkxz,Gxxzzr);


    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nzb,nxb);
    cur    =  sf_floatalloc2(nzb,nxb);
    aa     =  sf_floatalloc3(nzb,nxb,3);
    

    bd2_init(nx,nz,nxl,nxr,nbt,nbb,cxl,cxr,czt,czb);

	 lowrank_init2(nzb, nxb, nkxz, nkzz, nkxx, m2, nxzb2, Cxxl, Cxxr, Czzl, Czzr, Cxxxxl, Cxxxxr, Czzzzl, Czzzzr, Cxzzzl, Cxzzzr, Cxxxzl, Cxxxzr, Cxxzzl, Cxxzzr);

    /*input & extend velocity model*/
    v = sf_floatalloc2(nzb,nxb);
    vtmp = sf_floatalloc2(nz,nx);
    sf_floatread(vtmp[0],nx*nz,vel);
	v = extmodel(vtmp, nz, nx, nbt);

    sigma = sf_floatalloc2(nzb,nxb);
    sf_floatread(vtmp[0],nx*nz,fsigma);
	sigma = extmodel(vtmp, nz, nx, nbt);

    delta = sf_floatalloc2(nzb,nxb);
    sf_floatread(vtmp[0],nx*nz,fdelta);
	delta = extmodel(vtmp, nz, nx, nbt);
    
    seta = sf_floatalloc2(nzb,nxb);
    sf_floatread(vtmp[0],nx*nz,fseta);
	seta = extmodel(vtmp, nz, nx, nbt);

    v0 =0.0;
    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            v0 += v[ix][iz]*v[ix][iz];
         }
    }

    v0 = sqrtf(v0/(nxb*nzb));
	fprintf(stderr, "v0=%f\n\n", v0);

    for (ix=0; ix < nxb; ix++) {
         for (iz=0; iz < nzb; iz++) {
         w = v[ix][iz]*v[ix][iz];
         g1 = dt*dt*(v[ix][iz]*v[ix][iz]-v0*v0)/(12.0*dz*dz);
         g2 = dt*dt*(v[ix][iz]*v[ix][iz]-v0*v0)/(12.0*dx*dx);
         aa[1][ix][iz] = w*g1;
         aa[2][ix][iz] = w*g2;
         aa[0][ix][iz] = w-2.0*aa[1][ix][iz]-2.0*aa[2][ix][iz] ;
        }
    }

    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            cur[ix][iz] = 0.0;
            old[ix][iz] = 0.0; 
        }
    }

    /* propagation in time */
//	nkxz=fft2_init(true,1,nzb,nxb,&nkzz,&nkxx);
    for (it=0; it < nt; it++) {
		fprintf(stderr, "\b\b\b\b\b%d", it);

        opostep2(old, cur, nzb, nxb, dz, dx, v0, v, sigma, delta, seta, dt); 
        old[isx+nxl][isz+nbt] += wav[it];    
		bd2_decay(old); 
        bd2_decay(cur); 
        tmp = old;
        old = cur;
        cur = tmp;
        
		if (it%snap == 0)
		sf_floatwrite(cur[0], nxb*nzb, out);

    }
	lowrank_close2();
    bd2_close();

    free(**aa);
    free(*aa);
    free(aa);
    free(*v);     
    free(*sigma);     
    free(*delta);     
    free(*seta);     
	free(*vtmp);
    free(*cur);     
    free(*old);     
    free(v);     
    free(sigma);     
    free(delta);     
    free(seta);     
	free(vtmp);
    free(cur);     
    free(old);     
    exit(0); 
}           
           
