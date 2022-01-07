/* Modeling of pure acoustic wave in 3-D transversely isotropic meida using optimized pseudo-Laplacian operator */
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
#include "abcpass3.h"
#include "opostep3.h"
#include "fft3.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, ny, nz, nt, ix, iy, iz, it, nbt, nbb, nxl, nxr, nyl, nyr, nxb, nyb, nzb, isx, isy, isz;
    float dt, dx, dy, dz, o1, o2, o3;
    float ***old,  ***cur,  ***tmp, *wav;
    float  ***v, ***vtmp, v0, ***sigma, ***delta, ***seta, ***phi;
    float ***aa, czt, czb, cxl, cxr, cyl, cyr; /* top, bottom, left, right */
    float	**Cxxl, **Cxxr, **Cyyl, **Cyyr, **Czzl, **Czzr, **Cxyl, **Cxyr,
	**Cxzl, **Cxzr, **Cyzl, **Cyzr, **Cxxxxl, **Cxxxxr, **Cyyyyl,
	**Cyyyyr, **Czzzzl, **Czzzzr, **Cxxxyl, **Cxxxyr, **Cxxxzl,
	**Cxxxzr, **Cxyyyl, **Cxyyyr, **Cyyyzl, **Cyyyzr, **Cxzzzl,
	**Cxzzzr, **Cyzzzl, **Cyzzzr, **Cxxyyl, **Cxxyyr, **Cxxzzl,
	**Cxxzzr, **Cyyzzl, **Cyyzzr, **Cxxyzl, **Cxxyzr, **Cxyyzl,
	**Cxyyzr, **Cxyzzl, **Cxyzzr;
    sf_file out, vel, source, fsigma, fdelta, fseta, fphi,
	Gxxl, Gxxr, Gyyl, Gyyr, Gzzl, Gzzr, Gxyl, Gxyr,
	Gxzl, Gxzr, Gyzl, Gyzr, Gxxxxl, Gxxxxr, Gyyyyl,
	Gyyyyr, Gzzzzl, Gzzzzr, Gxxxyl, Gxxxyr, Gxxxzl,
	Gxxxzr, Gxyyyl, Gxyyyr, Gyyyzl, Gyyyzr, Gxzzzl,
	Gxzzzr, Gyzzzl, Gyzzzr, Gxxyyl, Gxxyyr, Gxxzzl,
	Gxxzzr, Gyyzzl, Gyyzzr, Gxxyzl, Gxxyzr, Gxyyzl,
	Gxyyzr, Gxyzzl, Gxyzzr;
    int opt, snap, nsnap;    /* optimal padding */
    int nkxyz,nkxx,nkyy,nkzz,nxyzb,nxyzb2,pad1,n2,m2;
    bool cmplx;
     
    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    fsigma = sf_input("sigma");   /* velocity */
    fdelta = sf_input("delta");   /* velocity */
    fseta  = sf_input("seta");   /* velocity */
    fphi  = sf_input("phi");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

    Gxxl  = sf_input("Gxxl");
    Gxxr  = sf_input("Gxxr");
    Gyyl  = sf_input("Gyyl");
    Gyyr  = sf_input("Gyyr");
    Gzzl  = sf_input("Gzzl");
    Gzzr  = sf_input("Gzzr");
	
    Gxyl  = sf_input("Gxyl");
    Gxyr  = sf_input("Gxyr");
    Gxzl  = sf_input("Gxzl");
    Gxzr  = sf_input("Gxzr");
    Gyzl  = sf_input("Gyzl");
    Gyzr  = sf_input("Gyzr");

    Gxxxxl  = sf_input("Gxxxxl");
    Gxxxxr  = sf_input("Gxxxxr");
    Gyyyyl  = sf_input("Gyyyyl");
    Gyyyyr  = sf_input("Gyyyyr");
    Gzzzzl  = sf_input("Gzzzzl");
    Gzzzzr  = sf_input("Gzzzzr");

    Gxxxyl  = sf_input("Gxxxyl");
    Gxxxyr  = sf_input("Gxxxyr");
    Gxxxzl  = sf_input("Gxxxzl");
    Gxxxzr  = sf_input("Gxxxzr");
    Gxyyyl  = sf_input("Gxyyyl");
    Gxyyyr  = sf_input("Gxyyyr");
    Gyyyzl  = sf_input("Gyyyzl");
    Gyyyzr  = sf_input("Gyyyzr");
    Gxzzzl  = sf_input("Gxzzzl");
    Gxzzzr  = sf_input("Gxzzzr");
    Gyzzzl  = sf_input("Gyzzzl");
    Gyzzzr  = sf_input("Gyzzzr");

    Gxxyyl  = sf_input("Gxxyyl");
    Gxxyyr  = sf_input("Gxxyyr");
    Gxxzzl  = sf_input("Gxxzzl");
    Gxxzzr  = sf_input("Gxxzzr");
    Gyyzzl  = sf_input("Gyyzzl");
    Gyyzzr  = sf_input("Gyyzzr");

    Gxxyzl  = sf_input("Gxxyzl");
    Gxxyzr  = sf_input("Gxxyzr");
    Gxyyzl  = sf_input("Gxyyzl");
    Gxyyzr  = sf_input("Gxyyzr");
    Gxyzzl  = sf_input("Gxyzzl");
    Gxyzzr  = sf_input("Gxyzzr");

    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");

    if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histint(vel,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histfloat(vel,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(vel,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(vel,"o2",&o2)) o2=0.0;
    if (!sf_histfloat(vel,"o3",&o3)) o3=0.0;

    if (!sf_getint("opt",&opt)) opt=1;
    /* if y, determine optimal size for efficiency */

    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isy",&isy)) sf_error("Need isy input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    if (!sf_getint("nbt",&nbt)) nbt=30;
    if (!sf_getint("nbb",&nbb)) nbb=30;
    if (!sf_getint("nxl",&nxl)) nxl=30;
    if (!sf_getint("nxr",&nxr)) nxr=30;
    if (!sf_getint("nyl",&nyl)) nyl=30;
    if (!sf_getint("nyr",&nyr)) nyr=30;
	
    /* assume ABC pars are the same */
    if (nbt != nbb || nxl != nxr || nyl != nyr || nbt!=nxl || nbt!=nyl) 
	sf_error("ABC pars are not the same");

    if (!sf_getfloat("czt",&czt))  czt = 0.01; /*decaying parameter*/
    if (!sf_getfloat("czb",&czb))  czb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxl",&cxl)) cxl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxr",&cxr)) cxr = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cyl",&cxl)) cyl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cyr",&cxr)) cyr = 0.01; /*decaying parameter*/

    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */

    if (!sf_getint("snap",&snap)) snap=1;
    nsnap=0;
    for (it=0; it < nt; it++) {
	if (it%snap == 0) nsnap++;
    }

    sf_putfloat(out,"d1",dz);
    sf_putfloat(out,"d2",dx);
    sf_putfloat(out,"d3",dy);
    sf_putfloat(out,"d4",dt*snap);
    sf_putfloat(out,"o1",o1); 
    sf_putfloat(out,"o2",o2); 
    sf_putfloat(out,"o3",o3); 
    sf_putfloat(out,"o4",0.0); 

    nxb = nx + nxl + nxr;
    nyb = ny + nyl + nyr;
    nzb = nz + nbt + nbb;
    
    sf_putint(out,"n1",nzb);
    sf_putint(out,"n2",nxb);
    sf_putint(out,"n3",nyb);
    sf_putint(out,"n4",nsnap);

    nkxyz=fft3_init(cmplx,pad1,nzb,nxb,nyb,&nkzz,&nkxx,&nkyy);
    nxyzb = nxb*nyb*nzb;
    nxyzb2 = nkxx*nkyy*nkzz;

    if (!sf_histint(Gxxl,"n1",&n2) || n2 != nxyzb) sf_error("Need n1=%d in left",nxyzb);
    if (!sf_histint(Gxxl,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(Gxxr,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(Gxxr,"n2",&n2) || n2 != nkxyz) sf_error("Need n2=%d in right",nkxyz);

    Cxxl = sf_floatalloc2(nxyzb,m2);
    Cxxr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxxl[0],nxyzb*m2,Gxxl);
    sf_floatread(Cxxr[0],m2*nkxyz,Gxxr);

    Cyyl = sf_floatalloc2(nxyzb,m2);
    Cyyr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cyyl[0],nxyzb*m2,Gyyl);
    sf_floatread(Cyyr[0],m2*nkxyz,Gyyr);
    
    Czzl = sf_floatalloc2(nxyzb,m2);
    Czzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Czzl[0],nxyzb*m2,Gzzl);
    sf_floatread(Czzr[0],m2*nkxyz,Gzzr);


    Cxyl = sf_floatalloc2(nxyzb,m2);
    Cxyr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxyl[0],nxyzb*m2,Gxyl);
    sf_floatread(Cxyr[0],m2*nkxyz,Gxyr);

    Cxzl = sf_floatalloc2(nxyzb,m2);
    Cxzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxzl[0],nxyzb*m2,Gxzl);
    sf_floatread(Cxzr[0],m2*nkxyz,Gxzr);

    Cyzl = sf_floatalloc2(nxyzb,m2);
    Cyzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cyzl[0],nxyzb*m2,Gyzl);
    sf_floatread(Cyzr[0],m2*nkxyz,Gyzr);


    Cxxxxl = sf_floatalloc2(nxyzb,m2);
    Cxxxxr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxxxxl[0],nxyzb*m2,Gxxxxl);
    sf_floatread(Cxxxxr[0],m2*nkxyz,Gxxxxr);

    Cyyyyl = sf_floatalloc2(nxyzb,m2);
    Cyyyyr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cyyyyl[0],nxyzb*m2,Gyyyyl);
    sf_floatread(Cyyyyr[0],m2*nkxyz,Gyyyyr);
    
    Czzzzl = sf_floatalloc2(nxyzb,m2);
    Czzzzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Czzzzl[0],nxyzb*m2,Gzzzzl);
    sf_floatread(Czzzzr[0],m2*nkxyz,Gzzzzr);

    Cxxxyl = sf_floatalloc2(nxyzb,m2);
    Cxxxyr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxxxyl[0],nxyzb*m2,Gxxxyl);
    sf_floatread(Cxxxyr[0],m2*nkxyz,Gxxxyr);

    Cxxxzl = sf_floatalloc2(nxyzb,m2);
    Cxxxzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxxxzl[0],nxyzb*m2,Gxxxzl);
    sf_floatread(Cxxxzr[0],m2*nkxyz,Gxxxzr);

    Cxyyyl = sf_floatalloc2(nxyzb,m2);
    Cxyyyr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxyyyl[0],nxyzb*m2,Gxyyyl);
    sf_floatread(Cxyyyr[0],m2*nkxyz,Gxyyyr);

    Cyyyzl = sf_floatalloc2(nxyzb,m2);
    Cyyyzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cyyyzl[0],nxyzb*m2,Gyyyzl);
    sf_floatread(Cyyyzr[0],m2*nkxyz,Gyyyzr);
	
    Cxzzzl = sf_floatalloc2(nxyzb,m2);
    Cxzzzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxzzzl[0],nxyzb*m2,Gxzzzl);
    sf_floatread(Cxzzzr[0],m2*nkxyz,Gxzzzr);

    Cyzzzl = sf_floatalloc2(nxyzb,m2);
    Cyzzzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cyzzzl[0],nxyzb*m2,Gyzzzl);
    sf_floatread(Cyzzzr[0],m2*nkxyz,Gyzzzr);

    Cxxyyl = sf_floatalloc2(nxyzb,m2);
    Cxxyyr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxxyyl[0],nxyzb*m2,Gxxyyl);
    sf_floatread(Cxxyyr[0],m2*nkxyz,Gxxyyr);

    Cxxzzl = sf_floatalloc2(nxyzb,m2);
    Cxxzzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxxzzl[0],nxyzb*m2,Gxxzzl);
    sf_floatread(Cxxzzr[0],m2*nkxyz,Gxxzzr);

    Cyyzzl = sf_floatalloc2(nxyzb,m2);
    Cyyzzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cyyzzl[0],nxyzb*m2,Gyyzzl);
    sf_floatread(Cyyzzr[0],m2*nkxyz,Gyyzzr);

    Cxxyzl = sf_floatalloc2(nxyzb,m2);
    Cxxyzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxxyzl[0],nxyzb*m2,Gxxyzl);
    sf_floatread(Cxxyzr[0],m2*nkxyz,Gxxyzr);

    Cxyyzl = sf_floatalloc2(nxyzb,m2);
    Cxyyzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxyyzl[0],nxyzb*m2,Gxyyzl);
    sf_floatread(Cxyyzr[0],m2*nkxyz,Gxyyzr);

    Cxyzzl = sf_floatalloc2(nxyzb,m2);
    Cxyzzr = sf_floatalloc2(m2,nkxyz);
    sf_floatread(Cxyzzl[0],nxyzb*m2,Gxyzzl);
    sf_floatread(Cxyzzr[0],m2*nkxyz,Gxyzzr);
   
   
    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc3(nzb,nxb,nyb);
    cur    =  sf_floatalloc3(nzb,nxb,nyb);
    aa     =  sf_floatalloc3(nzb,nxb,3);
    

    bd3_init(ny,nx,nz,nyl,nyr,nxl,nxr,nbt,nbb,cyl,cyr,cxl,cxr,czt,czb);
    lowrank_init3(nzb, nxb, nyb, nkxyz, nkzz, nkxx, nkyy, m2, nxyzb2, Cxxl, Cxxr, Cyyl, Cyyr, Czzl, Czzr, Cxyl, Cxyr, Cxzl, Cxzr, Cyzl, Cyzr, Cxxxxl, 
		  Cxxxxr, Cyyyyl, Cyyyyr, Czzzzl, Czzzzr, Cxxxyl, Cxxxyr, Cxxxzl, Cxxxzr, Cxyyyl, Cxyyyr, Cyyyzl, Cyyyzr, Cxzzzl, Cxzzzr, Cyzzzl,
		  Cyzzzr, Cxxyyl, Cxxyyr, Cxxzzl, Cxxzzr, Cyyzzl, Cyyzzr, Cxxyzl, Cxxyzr, Cxyyzl, Cxyyzr, Cxyzzl, Cxyzzr);

    /*input & extend velocity model*/
    v = sf_floatalloc3(nzb,nxb,nyb);
    vtmp = sf_floatalloc3(nz,nx,ny);
    sf_floatread(vtmp[0][0],nx*ny*nz,vel);
    v = extmodel3d(vtmp, nz, nx, ny, nbt);

    sigma = sf_floatalloc3(nzb,nxb,nyb);
    sf_floatread(vtmp[0][0],nx*ny*nz,fsigma);
    sigma = extmodel3d(vtmp, nz, nx, ny, nbt);

    delta = sf_floatalloc3(nzb,nxb,nyb);
    sf_floatread(vtmp[0][0],nx*ny*nz,fdelta);
    delta = extmodel3d(vtmp, nz, nx, ny, nbt);
    
    seta = sf_floatalloc3(nzb,nxb,nyb);
    sf_floatread(vtmp[0][0],nx*ny*nz,fseta);
    seta = extmodel3d(vtmp, nz, nx, ny, nbt);
    
    phi = sf_floatalloc3(nzb,nxb,nyb);
    sf_floatread(vtmp[0][0],nx*ny*nz,fphi);
    phi = extmodel3d(vtmp, nz, nx, ny, nbt);

    v0 =0.0;
    for (iy=0; iy < nyb; iy++) {
	for (ix=0; ix < nxb; ix++) {
	    for (iz=0; iz < nzb; iz++) {
		v0 += v[iy][ix][iz]*v[iy][ix][iz];
	    }
	}
    }

    v0 = sqrtf(v0/(nxb*nyb*nzb));

    fprintf(stderr, "v0=%f\n\n", v0);


    for (iy=0; iy < nyb; iy++) {
	for (ix=0; ix < nxb; ix++) {
	    for (iz=0; iz < nzb; iz++) {
		cur[iy][ix][iz] = 0.0;
		old[iy][ix][iz] = 0.0; 
	    }
	}
    }

    /* propagation in time */
//	nkxyz=fft3_init(true,1,nzb,nxb,nyb,&nkzz,&nkxx,&nkyy);
    fprintf(stderr, "nkzz=%d nkxx=%d nkyy=%d\n", nkzz, nkxx, nkyy);
    for (it=0; it < nt; it++) {
	fprintf(stderr, "\b\b\b\b\b%d", it);

        opostep3(old, cur, nzb, nxb, nyb, dz, dx, dy, v0, v, sigma, delta, seta, phi, dt); 
        old[isy+nyl][isx+nxl][isz+nbt] += wav[it];
      
	bd3_decay(old); 
        bd3_decay(cur); 
        tmp = old;
        old = cur;
        cur = tmp;
	if (it%snap == 0)
	    sf_floatwrite(cur[0][0], nxb*nyb*nzb, out);

    }
    lowrank_close3();
    bd3_close();

    free(**aa);
    free(*aa);
    free(aa);
    free(**v);     
    free(**sigma);     
    free(**delta);     
    free(**seta);     
    free(**vtmp);
    free(**cur);     
    free(**old);     
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
