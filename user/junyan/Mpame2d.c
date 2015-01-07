/*  2-D elasitc wave modeling and vector field decompostion using pseudo-analytical method */
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
#include "abcpass2.h"
#include "pamstep2e.h"
#include "fft2.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nz, nt, ix, iz, it, nbt, nbb, nxl, nxr,  nxb, nyb, nzb, isx, isz;
    float dt, dx, dy, dz, o1, o2, o3;
    float **upold,  **upcur, **wpold, **wpcur, *wav;
	 float **usold, **uscur, **wsold, **wscur, **tmp;
	 float **uu,**ww;
    float  **vp, vp0, **vs, vs0, **vtmp ;
    float ***aa, w, g1, g2, czt, czb, cxl, cxr; /* top, bottom, left, right */
    float ax, az, factor;
    sf_file fwavup, fwavwp, fwavus, fwavws, fwavu, fwavw, fvelp, fvels, fsource ;
    int opt, snap, nsnap;    /* optimal padding */
	 int nkxz,nkxx,nkzz;
     
    sf_init(argc,argv);
    fvelp = sf_input("vp");   /* velocity */
    fvels = sf_input("vs");   /* velocity */
    fsource = sf_input("in");   /* source wavlet*/
    fwavup = sf_output("out");
    fwavwp = sf_output("wavwp");
    fwavus = sf_output("wavus");
    fwavws = sf_output("wavws");
    fwavu = sf_output("wavu");
    fwavw = sf_output("wavw");

    if (SF_FLOAT != sf_gettype(fvelp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fvels)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fsource)) sf_error("Need float input");

	if (!sf_histint(fvelp,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(fvelp,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(fvelp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(fvelp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(fvelp,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(fvelp,"o2",&o2)) o2=0.0;

    
	if (!sf_getint("opt",&opt)) opt=0;
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

    sf_putfloat(fwavup,"d1",dz);
    sf_putfloat(fwavup,"d2",dx);
    sf_putfloat(fwavup,"d3",dt*snap);
    sf_putfloat(fwavup,"o1",o1); 
    sf_putfloat(fwavup,"o2",o2); 
    sf_putfloat(fwavup,"o3",0.0);

    sf_putfloat(fwavwp,"d1",dz);
    sf_putfloat(fwavwp,"d2",dx);
    sf_putfloat(fwavwp,"d3",dt*snap);
    sf_putfloat(fwavwp,"o1",o1); 
    sf_putfloat(fwavwp,"o2",o2); 
    sf_putfloat(fwavwp,"o3",0.0);

    sf_putfloat(fwavus,"d1",dz);
    sf_putfloat(fwavus,"d2",dx);
    sf_putfloat(fwavus,"d3",dt*snap);
    sf_putfloat(fwavus,"o1",o1); 
    sf_putfloat(fwavus,"o2",o2); 
    sf_putfloat(fwavus,"o3",0.0);

    sf_putfloat(fwavws,"d1",dz);
    sf_putfloat(fwavws,"d2",dx);
    sf_putfloat(fwavws,"d3",dt*snap);
    sf_putfloat(fwavws,"o1",o1); 
    sf_putfloat(fwavws,"o2",o2); 
    sf_putfloat(fwavws,"o3",0.0);


    sf_putfloat(fwavu,"d1",dz);
    sf_putfloat(fwavu,"d2",dx);
    sf_putfloat(fwavu,"d3",dt*snap);
    sf_putfloat(fwavu,"o1",o1); 
    sf_putfloat(fwavu,"o2",o2); 
    sf_putfloat(fwavu,"o3",0.0);

    sf_putfloat(fwavw,"d1",dz);
    sf_putfloat(fwavw,"d2",dx);
    sf_putfloat(fwavw,"d3",dt*snap);
    sf_putfloat(fwavw,"o1",o1); 
    sf_putfloat(fwavw,"o2",o2); 
    sf_putfloat(fwavw,"o3",0.0);

    nxb = nx + nxl + nxr;
    nzb = nz + nbt + nbb;
    
	sf_putint(fwavup,"n1",nzb);
    sf_putint(fwavup,"n2",nxb);
    sf_putint(fwavup,"n3",nsnap);

	sf_putint(fwavwp,"n1",nzb);
    sf_putint(fwavwp,"n2",nxb);
    sf_putint(fwavwp,"n3",nsnap);

	sf_putint(fwavus,"n1",nzb);
    sf_putint(fwavus,"n2",nxb);
    sf_putint(fwavus,"n3",nsnap);

	sf_putint(fwavws,"n1",nzb);
    sf_putint(fwavws,"n2",nxb);
    sf_putint(fwavws,"n3",nsnap);

	sf_putint(fwavu,"n1",nzb);
    sf_putint(fwavu,"n2",nxb);
    sf_putint(fwavu,"n3",nsnap);


	sf_putint(fwavw,"n1",nzb);
    sf_putint(fwavw,"n2",nxb);
    sf_putint(fwavw,"n3",nsnap);



    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,fsource);

    upold    =  sf_floatalloc2(nzb,nxb);
    upcur    =  sf_floatalloc2(nzb,nxb);
    wpold    =  sf_floatalloc2(nzb,nxb);
    wpcur    =  sf_floatalloc2(nzb,nxb);

    usold    =  sf_floatalloc2(nzb,nxb);
    uscur    =  sf_floatalloc2(nzb,nxb);
    wsold    =  sf_floatalloc2(nzb,nxb);
    wscur    =  sf_floatalloc2(nzb,nxb);

    uu   =  sf_floatalloc2(nzb,nxb);
    ww    =  sf_floatalloc2(nzb,nxb);

    aa     =  sf_floatalloc3(nzb,nxb,3);
    

    bd2_init(nx,nz,nxl,nxr,nbt,nbb,cxl,cxr,czt,czb);

    /*input & extend velocity model*/
    vp = sf_floatalloc2(nzb,nxb);
    vtmp = sf_floatalloc2(nz,nx);
    sf_floatread(vtmp[0],nx*nz,fvelp);
	 vp = extmodel(vtmp, nz, nx, nbt);

    vs = sf_floatalloc2(nzb,nxb);
    sf_floatread(vtmp[0],nx*nz,fvels);
	 vs= extmodel(vtmp, nz, nx, nbt);


    vp0 =0.0;
    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            vp0 += vp[ix][iz]*vp[ix][iz];
         }
    }

    vp0 = sqrtf(vp0/(nxb*nzb));

	fprintf(stderr, "vp0=%f\n\n", vp0);

   vs0 =0.0;
    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            vs0 += vs[ix][iz]*vs[ix][iz];
         }
    }

    vs0 = sqrtf(vs0/(nxb*nzb));

	fprintf(stderr, "vs0=%f\n\n", vs0);


    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            upcur[ix][iz] = 0.0;
            upold[ix][iz] = 0.0; 
            wpcur[ix][iz] = 0.0;
            wpold[ix][iz] = 0.0; 
            uscur[ix][iz] = 0.0;
            usold[ix][iz] = 0.0; 
            wscur[ix][iz] = 0.0;
            wsold[ix][iz] = 0.0; 
				
			uu[ix][iz] = 0.0;
			ww[ix][iz] = 0.0;
        }
    }

    /* propagation in time */
    pamstep2e_init(nzb,nxb,dz,dx,opt);

	nkxz=fft2_init(true,1,nzb,nxb,&nkzz,&nkxx);
    for (it=0; it < nt; it++) {
		fprintf(stderr, "\b\b\b\b\b%d", it);

   //     uu[isx+nxl][isz+nbt] += wav[it];
   //     ww[isx+nxl][isz+nbt] += wav[it];

        pamstep2e(upold, upcur, wpold, wpcur, usold, uscur, wsold, wscur, uu, ww, nzb, nxb, dz, dx, vp0, vs0, vp, vs, dt); 
   	 
	 for (ix=0; ix < nxb; ix++) {
         for (iz=0; iz < nzb; iz++) {
				uu[ix][iz] = upold[ix][iz] + usold[ix][iz];
				ww[ix][iz] = wpold[ix][iz] + wsold[ix][iz];
		}
	}
			
      
   	 uu[isx+nxl][isz+nbt] += wav[it];
//        ww[isx+nxl][isz+nbt] += wav[it];

		bd2_decay(upold); 
		bd2_decay(upcur); 
		bd2_decay(usold); 
        bd2_decay(uscur); 
		bd2_decay(wpold); 
        bd2_decay(wpcur); 
		bd2_decay(wsold); 
        bd2_decay(wscur); 


        tmp = upold;
        upold = upcur;
        upcur = tmp;

        tmp = wpold;
        wpold = wpcur;
        wpcur = tmp;

        tmp = usold;
        usold = uscur;
        uscur = tmp;

        tmp = wsold;
        wsold = wscur;
        wscur = tmp;

        
       /* for (iz=nbt; iz<nz+nbt; iz++){
             sf_floatwrite(cur[iz]+nbl,nx,out);
         }
	   */
		if (it%snap==0) {
		sf_floatwrite(upcur[0], nxb*nzb, fwavup);
		sf_floatwrite(wpcur[0], nxb*nzb, fwavwp);
		sf_floatwrite(uscur[0], nxb*nzb, fwavus);
		sf_floatwrite(wscur[0], nxb*nzb, fwavws);
		sf_floatwrite(uu[0], nxb*nzb, fwavu);
		sf_floatwrite(ww[0], nxb*nzb, fwavw);
		}

    }
    pamstep2e_close();
    bd2_close();

    free(**aa);
    free(*aa);
    free(aa);
    free(*vp);     
     
	free(*vtmp);
    free(*upcur);     
    free(*upold);     
    free(vp);     
    
	free(vtmp);
    free(upcur);     
    free(upold);     
    exit(0); 
}           
           
