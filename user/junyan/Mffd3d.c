/* 3-D Fourier finite-difference wave extrapolation */
/*
  Copyright (C) 2014 Institute of Geology and Geophysics, Chinese Academy of Sciences (Jun Yan) 
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
#include "ffdstep3.h"
#include "srcsm.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, ny, nz, nt, ix, iy, iz, it, nbt, nbb, nxl, nxr, nyl, nyr, nxb, nyb, nzb, isx, isy, isz;
    float dt, dx, dy, dz, o1, o2, o3;
    float ***old,  ***cur,  ***tmp, *wav;
    float  ***v, ***vtmp, v0, ****aa, w, g1, g2, g3, czt, czb, cxl, cxr, cyl, cyr; /* top, bottom, left, right */
    float ax, ay, az, factor;
    sf_file out, vel, source;
    bool opt;    /* optimal padding */
     
    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
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
    /*  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input"); */
    /*  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input"); */
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */

    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isy",&isy)) sf_error("Need isy input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    if (!sf_getint("nbt",&nbt)) nbt=44;
    if (!sf_getint("nbb",&nbb)) nbb=44;
    if (!sf_getint("nxl",&nxl)) nxl=44;
    if (!sf_getint("nxr",&nxr)) nxr=44;
    if (!sf_getint("nyl",&nyl)) nyl=44;
    if (!sf_getint("nyr",&nyr)) nyr=44;
	
	/* assume ABC pars are the same */
	if (nbt != nbb || nxl != nxr || nyl != nyr || nbt!=nxl || nbt!=nyl) 
		sf_error("ABC pars are not the same");

    if (!sf_getfloat("czt",&czt))  czt = 0.01; /*decaying parameter*/
    if (!sf_getfloat("czb",&czb))  czb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxl",&cxl)) cxl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxr",&cxr)) cxr = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cyl",&cyl)) cyl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cyr",&cyr)) cyr = 0.01; /*decaying parameter*/

    if (!sf_getfloat("ax",&ax)) ax= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("ay",&ay)) ay= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("az",&az)) az= 5.0; /*suppress HF parameter*/
    if (!sf_getfloat("factor",&factor)) factor= 5.0/6.0; /*suppress HF parameter*/

    sf_putfloat(out,"d1",dz);
    sf_putfloat(out,"d2",dx);
    sf_putfloat(out,"d3",dy);
    sf_putfloat(out,"d4",dt);
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
    sf_putint(out,"n4",nt);


    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc3(nzb,nxb,nyb);
    cur    =  sf_floatalloc3(nzb,nxb,nyb);
    aa     =  sf_floatalloc4(nzb,nxb,nyb,4);
    

    bd3_init(ny,nx,nz,nyl,nyr,nxl,nxr,nbt,nbb,cyl,cyr,cxl,cxr,czt,czb);

    /*input & extend velocity model*/
    v = sf_floatalloc3(nzb,nxb,nyb);
    vtmp = sf_floatalloc3(nz,nx,ny);
	
    sf_floatread(vtmp[0][0],nx*ny*nz,vel);


	v = extmodel3d(vtmp, nz, nx, ny, nbt);

    v0 =0.0;
    for (iy=0; iy < nyb; iy++) {
        for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            v0 += v[iy][ix][iz]*v[iy][ix][iz];
         }
    }
	}

    v0 = sqrtf(v0/(nxb*nyb*nzb));

    for (iy=0; iy < nyb; iy++){
         for (ix=0; ix < nxb; ix++) {
         for (iz=0; iz < nzb; iz++) {
         w = v[iy][ix][iz]*v[iy][ix][iz];
         g1 = dt*dt*(v[iy][ix][iz]*v[iy][ix][iz]-v0*v0)/(12.0*dz*dz);
         g2 = dt*dt*(v[iy][ix][iz]*v[iy][ix][iz]-v0*v0)/(12.0*dx*dx);
         g3 = dt*dt*(v[iy][ix][iz]*v[iy][ix][iz]-v0*v0)/(12.0*dy*dy);
         aa[1][iy][ix][iz] = w*g1;
         aa[2][iy][ix][iz] = w*g2;
         aa[3][iy][ix][iz] = w*g3; 
         aa[0][iy][ix][iz] = w-2.0*aa[1][iy][ix][iz]-2.0*aa[2][iy][ix][iz]-2.0*aa[3][iy][ix][iz] ;
        }
      }
	}

    for (iy=0; iy < nyb; iy++) {
        for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            cur[iy][ix][iz] = 0.0;
            old[iy][ix][iz] = 0.0; 
        }
    }
	}

    /* propagation in time */
    ffdstep3_init(nzb,nxb,nyb,dz,dx,dy);
//  srcsm_init(dz,dx);

    for (it=0; it < nt; it++) {
		fprintf(stderr, "\b\b\b\b\b%d", it);

        ffdstep3_dehf(old, cur, aa, nzb, nxb, nyb, v0, dt, az, ax, ay, factor); 
        old[isy+nyl][isx+nxl][isz+nbt] += wav[it];
//      source_smooth(old,isz+nbt,isx+nbl,wav[it]);
        bd3_decay(old); 
        bd3_decay(cur); 
        tmp = old;
        old = cur;
        cur = tmp;
        
       /* for (iz=nbt; iz<nz+nbt; iz++){
             sf_floatwrite(cur[iz]+nbl,nx,out);
         }
	   */
		
		sf_floatwrite(cur[0][0], nxb*nyb*nzb, out);

    }

    ffdstep3_close();
    bd3_close();

    free(***aa);
    free(**aa);
    free(*aa);
    free(aa);
    free(**v);
	free(**vtmp);
    free(**cur);     
    free(**old);     
    free(*v);     
	free(*vtmp);
    free(*cur);     
    free(*old);     
    free(v);     
	free(vtmp);
    free(cur);     
    free(old);     
    exit(0); 
}           
           
