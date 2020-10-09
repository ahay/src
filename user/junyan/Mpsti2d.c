/* Modeling of pure acoustic wave in 2-D transversely isotropic meida using psuedospectral method */
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
#include "psstep2.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nz, nt, ix, iz, it, nbt, nbb, nxl, nxr,  nxb, nzb, isx, isz;
    float dt, dx, dy, dz, o1, o2, o3;
    float **old,  **cur,  **tmp, *wav;
    float  **v, **vtmp, v0, **sigma, **delta, **seta;
    float ***aa, w, g1, g2, czt, czb, cxl, cxr; /* top, bottom, left, right */
    float ax, az, factor;
    sf_file out, vel, source, fsigma, fdelta, fseta;
    int opt, snap, nsnap;    /* optimal padding */
     
    sf_init(argc,argv);
    out = sf_output("out");
    vel = sf_input("vel");   /* velocity */
    fsigma = sf_input("sigma");   /* velocity */
    fdelta = sf_input("delta");   /* velocity */
    fseta  = sf_input("seta");   /* velocity */
    source = sf_input("in");   /* source wavlet*/

    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");

	if (!sf_histint(vel,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(vel,"o2",&o2)) o2=0.0;

    
	if (!sf_getint("opt",&opt)) opt=1;
    /* if y, determine optimal size for efficiency */

    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isx",&isx)) sf_error("Need isx input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");

    if (!sf_getint("nbt",&nbt)) nbt=44;
    if (!sf_getint("nbb",&nbb)) nbb=44;
    if (!sf_getint("nxl",&nxl)) nxl=44;
    if (!sf_getint("nxr",&nxr)) nxr=44;
	
	/* assume ABC pars are the same */
	if (nbt != nbb || nxl != nxr || nbt!=nxl) 
		sf_error("ABC pars are not the same");

    if (!sf_getfloat("czt",&czt))  czt = 0.01; /*decaying parameter*/
    if (!sf_getfloat("czb",&czb))  czb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxl",&cxl)) cxl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxr",&cxr)) cxr = 0.01; /*decaying parameter*/
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


    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,source);

    old    =  sf_floatalloc2(nzb,nxb);
    cur    =  sf_floatalloc2(nzb,nxb);
    aa     =  sf_floatalloc3(nzb,nxb,3);
    

    bd2_init(nx,nz,nxl,nxr,nbt,nbb,cxl,cxr,czt,czb);

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
            cur[ix][iz] = 0.0;
            old[ix][iz] = 0.0; 
        }
    }

    /* propagation in time */
    psstep2_init(nzb,nxb,dz,dx,opt);

    for (it=0; it < nt; it++) {
		fprintf(stderr, "\b\b\b\b\b%d", it);

        psstep2(old, cur, nzb, nxb, dz, dx, v0, v, sigma, delta, seta, dt); 
        old[isx+nxl][isz+nbt] += wav[it];
    
		  bd2_decay(old); 
        bd2_decay(cur); 
        tmp = old;
        old = cur;
        cur = tmp;
        

		if (it%snap==0)
		sf_floatwrite(cur[0], nxb*nzb, out);

    }

    psstep2_close();
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
           
