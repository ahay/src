/*  3-D elasitc wave modeling and vector field decompostion using pseudospectra method */
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
#include "psstep3e.h"
//#include "srcsm2.h"
#include "fft3.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, ny, nz, nt, ix, iy, iz, it, nbt, nbb, nxl, nxr, nyl, nyr,  nxb, nyb, nzb, isx, isy, isz;
    float dt, dx, dy, dz, o1, o2, o3, o4;
    float ***upold,  ***upcur, ***vpold, ***vpcur,  ***wpold, ***wpcur, *wav;
	 float ***usold, ***uscur, ***vsold, ***vscur, ***wsold, ***wscur, ***tmp;
	 float ***uu,***vv,***ww;
    float  ***vp, vp0, ***vs, vs0, ***vtmp ;
    float ***aa, w, g1, g2, g3, czt, czb, cxl, cxr, cyl, cyr; /* top, bottom, left, right */
    float ax, ay, az, factor;
    sf_file fwavup, fwavvp, fwavwp, fwavus, fwavvs, fwavws, fwavu, fwavv, fwavw, fvelp, fvels, fsource ;
    int opt, snap, nsnap;    /* optimal padding */
	 int nkxyz,nkxx,nkyy,nkzz;
     
    sf_init(argc,argv);
    fvelp = sf_input("vp");   /* velocity */
    fvels = sf_input("vs");   /* velocity */
    fsource = sf_input("in");   /* source wavlet*/
    fwavup = sf_output("out");
    fwavvp = sf_output("wavvp");
    fwavwp = sf_output("wavwp");
    fwavus = sf_output("wavus");
    fwavvs = sf_output("wavvs");
    fwavws = sf_output("wavws");
    fwavu = sf_output("wavu");
    fwavv = sf_output("wavv");
    fwavw = sf_output("wavw");

/*    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input"); */
    if (SF_FLOAT != sf_gettype(fvelp)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fvels)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(fsource)) sf_error("Need float input");

	if (!sf_histint(fvelp,"n1",&nz)) sf_error("No n1= in input");
	if (!sf_histfloat(fvelp,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(fvelp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(fvelp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histint(fvelp,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histfloat(fvelp,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(fvelp,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(fvelp,"o2",&o2)) o2=0.0;
    if (!sf_histfloat(fvelp,"o3",&o3)) o3=0.0;
    /*  if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input"); */
    /*  if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2= in input"); */
    
	if (!sf_getint("opt",&opt)) opt=1;
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
	if (nbt != nbb || nxl != nxr ||nyl!=nyr|| nbt!=nxl || nbt!=nyl) 
		sf_error("ABC pars are not the same");

    if (!sf_getfloat("czt",&czt))  czt = 0.01; /*decaying parameter*/
    if (!sf_getfloat("czb",&czb))  czb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxl",&cxl)) cxl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cxr",&cxr)) cxr = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cyl",&cxl)) cyl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cyr",&cxr)) cyr = 0.01; /*decaying parameter*/

	 if (!sf_getint("snap",&snap)) snap=1;
	 nsnap=0;
	 for (it=0; it < nt; it++) {
			if (it%snap == 0) nsnap++;
	  }

    sf_putfloat(fwavup,"d1",dz);
    sf_putfloat(fwavup,"d2",dx);
    sf_putfloat(fwavup,"d3",dy);
    sf_putfloat(fwavup,"d4",dt*snap);
    sf_putfloat(fwavup,"o1",o1); 
    sf_putfloat(fwavup,"o2",o2); 
    sf_putfloat(fwavup,"o3",o3); 
    sf_putfloat(fwavup,"o4",0.0);

    sf_putfloat(fwavvp,"d1",dz);
    sf_putfloat(fwavvp,"d2",dx);
    sf_putfloat(fwavvp,"d3",dy);
    sf_putfloat(fwavvp,"d4",dt*snap);
    sf_putfloat(fwavvp,"o1",o1); 
    sf_putfloat(fwavvp,"o2",o2); 
    sf_putfloat(fwavvp,"o3",o3); 
    sf_putfloat(fwavvp,"o4",0.0);
    
	sf_putfloat(fwavwp,"d1",dz);
    sf_putfloat(fwavwp,"d2",dx);
    sf_putfloat(fwavwp,"d3",dy);
    sf_putfloat(fwavwp,"d4",dt*snap);
    sf_putfloat(fwavwp,"o1",o1); 
    sf_putfloat(fwavwp,"o2",o2); 
    sf_putfloat(fwavwp,"o3",o3); 
    sf_putfloat(fwavwp,"o4",0.0);

    sf_putfloat(fwavus,"d1",dz);
    sf_putfloat(fwavus,"d2",dx);
    sf_putfloat(fwavus,"d3",dy);
    sf_putfloat(fwavus,"d4",dt*snap);
    sf_putfloat(fwavus,"o1",o1); 
    sf_putfloat(fwavus,"o2",o2); 
    sf_putfloat(fwavus,"o3",o3); 
    sf_putfloat(fwavus,"o4",0.0);

    sf_putfloat(fwavvs,"d1",dz);
    sf_putfloat(fwavvs,"d2",dx);
    sf_putfloat(fwavvs,"d3",dy);
    sf_putfloat(fwavvs,"d4",dt*snap);
    sf_putfloat(fwavvs,"o1",o1); 
    sf_putfloat(fwavvs,"o2",o2); 
    sf_putfloat(fwavvs,"o3",o3); 
    sf_putfloat(fwavvs,"o4",0.0);

    sf_putfloat(fwavws,"d1",dz);
    sf_putfloat(fwavws,"d2",dx);
    sf_putfloat(fwavws,"d3",dy);
    sf_putfloat(fwavws,"d4",dt*snap);
    sf_putfloat(fwavws,"o1",o1); 
    sf_putfloat(fwavws,"o2",o2); 
    sf_putfloat(fwavws,"o3",o3); 
    sf_putfloat(fwavws,"o4",0.0);


    sf_putfloat(fwavu,"d1",dz);
    sf_putfloat(fwavu,"d2",dx);
    sf_putfloat(fwavu,"d3",dy);
    sf_putfloat(fwavu,"d4",dt*snap);
    sf_putfloat(fwavu,"o1",o1); 
    sf_putfloat(fwavu,"o2",o2); 
    sf_putfloat(fwavu,"o3",o3); 
    sf_putfloat(fwavu,"o4",0.0);


    sf_putfloat(fwavv,"d1",dz);
    sf_putfloat(fwavv,"d2",dx);
    sf_putfloat(fwavv,"d3",dy);
    sf_putfloat(fwavv,"d4",dt*snap);
    sf_putfloat(fwavv,"o1",o1); 
    sf_putfloat(fwavv,"o2",o2); 
    sf_putfloat(fwavv,"o3",o3); 
    sf_putfloat(fwavv,"o4",0.0);

    sf_putfloat(fwavw,"d1",dz);
    sf_putfloat(fwavw,"d2",dx);
    sf_putfloat(fwavw,"d3",dy);
    sf_putfloat(fwavw,"d4",dt*snap);
    sf_putfloat(fwavw,"o1",o1); 
    sf_putfloat(fwavw,"o2",o2); 
    sf_putfloat(fwavw,"o3",o3); 
    sf_putfloat(fwavw,"o4",0.0);

    nxb = nx + nxl + nxr;
    nyb = ny + nyl + nyr;
    nzb = nz + nbt + nbb;
    
	sf_putint(fwavup,"n1",nzb);
    sf_putint(fwavup,"n2",nxb);
    sf_putint(fwavup,"n3",nyb);
    sf_putint(fwavup,"n4",nsnap);

	sf_putint(fwavvp,"n1",nzb);
    sf_putint(fwavvp,"n2",nxb);
    sf_putint(fwavvp,"n3",nyb);
    sf_putint(fwavvp,"n4",nsnap);

	sf_putint(fwavwp,"n1",nzb);
    sf_putint(fwavwp,"n2",nxb);
    sf_putint(fwavwp,"n3",nyb);
    sf_putint(fwavwp,"n4",nsnap);

	sf_putint(fwavus,"n1",nzb);
    sf_putint(fwavus,"n2",nxb);
    sf_putint(fwavus,"n3",nyb);
    sf_putint(fwavus,"n4",nsnap);

	sf_putint(fwavvs,"n1",nzb);
    sf_putint(fwavvs,"n2",nxb);
    sf_putint(fwavvs,"n3",nyb);
    sf_putint(fwavvs,"n4",nsnap);

	sf_putint(fwavws,"n1",nzb);
    sf_putint(fwavws,"n2",nxb);
    sf_putint(fwavws,"n3",nyb);
    sf_putint(fwavws,"n4",nsnap);

	sf_putint(fwavu,"n1",nzb);
    sf_putint(fwavu,"n2",nxb);
    sf_putint(fwavu,"n3",nyb);
    sf_putint(fwavu,"n4",nsnap);

	sf_putint(fwavv,"n1",nzb);
    sf_putint(fwavv,"n2",nxb);
    sf_putint(fwavv,"n3",nyb);
    sf_putint(fwavv,"n4",nsnap);

	sf_putint(fwavw,"n1",nzb);
    sf_putint(fwavw,"n2",nxb);
    sf_putint(fwavw,"n3",nyb);
    sf_putint(fwavw,"n4",nsnap);


    wav    =  sf_floatalloc(nt);
    sf_floatread(wav,nt,fsource);

    upold    =  sf_floatalloc3(nzb,nxb,nyb);
    upcur    =  sf_floatalloc3(nzb,nxb,nyb);
    vpold    =  sf_floatalloc3(nzb,nxb,nyb);
    vpcur    =  sf_floatalloc3(nzb,nxb,nyb);
    wpold    =  sf_floatalloc3(nzb,nxb,nyb);
    wpcur    =  sf_floatalloc3(nzb,nxb,nyb);

    usold    =  sf_floatalloc3(nzb,nxb,nyb);
    uscur    =  sf_floatalloc3(nzb,nxb,nyb);
    vsold    =  sf_floatalloc3(nzb,nxb,nyb);
    vscur    =  sf_floatalloc3(nzb,nxb,nyb);
    wsold    =  sf_floatalloc3(nzb,nxb,nyb);
    wscur    =  sf_floatalloc3(nzb,nxb,nyb);

    uu   =  sf_floatalloc3(nzb,nxb,nyb);
    vv   =  sf_floatalloc3(nzb,nxb,nyb);
    ww    =  sf_floatalloc3(nzb,nxb,nyb);

    aa     =  sf_floatalloc3(nzb,nxb,3);
    

    bd3_init(ny,nx,nz,nyl,nyr,nxl,nxr,nbt,nbb,cyl,cyr,cxl,cxr,czt,czb);

    /*input & extend velocity model*/
    vp = sf_floatalloc3(nzb,nxb,nyb);
    vtmp = sf_floatalloc3(nz,nx,ny);
    sf_floatread(vtmp[0][0],nx*ny*nz,fvelp);
	 vp = extmodel3d(vtmp, nz, nx, ny, nbt);

    vs = sf_floatalloc3(nzb,nxb,nyb);
    sf_floatread(vtmp[0][0],nx*ny*nz,fvels);
	 vs= extmodel3d(vtmp, nz, nx, ny, nbt);


    vp0 =0.0;
    for (iy=0; iy < nyb; iy++) {
    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            vp0 += vp[iy][ix][iz]*vp[iy][ix][iz];
         }
    }
	}

    vp0 = sqrtf(vp0/(nxb*nyb*nzb));
	fprintf(stderr, "vp0=%f\n\n", vp0);

	vs0 =0.0;
    for (iy=0; iy < nyb; iy++) {
    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            vs0 += vs[iy][ix][iz]*vs[iy][ix][iz];
         }
    }
	}

    vs0 = sqrtf(vs0/(nxb*nyb*nzb));
	fprintf(stderr, "vs0=%f\n\n", vs0);


    for (iy=0; iy < nyb; iy++) {
    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            upcur[iy][ix][iz] = 0.0;
            upold[iy][ix][iz] = 0.0; 
            vpcur[iy][ix][iz] = 0.0;
            vpold[iy][ix][iz] = 0.0; 
            wpcur[iy][ix][iz] = 0.0;
            wpold[iy][ix][iz] = 0.0; 
            uscur[iy][ix][iz] = 0.0;
            usold[iy][ix][iz] = 0.0; 
            vscur[iy][ix][iz] = 0.0;
            vsold[iy][ix][iz] = 0.0; 
            wscur[iy][ix][iz] = 0.0;
            wsold[iy][ix][iz] = 0.0; 
				
			uu[iy][ix][iz] = 0.0;
			vv[iy][ix][iz] = 0.0;
			ww[iy][ix][iz] = 0.0;
        }
    }
	}

    /* propagation in time */
    psstep3e_init(nzb,nxb,nyb,dz,dx,dy,opt);

	nkxyz=fft3_init(true,1,nzb,nxb,nyb,&nkzz,&nkxx,&nkyy);
    for (it=0; it < nt; it++) {
		fprintf(stderr, "\b\b\b\b\b%d", it);

   //     uu[isx+nxl][isz+nbt] += wav[it];
   //     ww[isx+nxl][isz+nbt] += wav[it];

        psstep3e(upold, upcur, vpold, vpcur, wpold, wpcur, usold, uscur, vsold, vscur, wsold, wscur, uu, vv, ww, nzb, nxb, nyb, dz, dx, dy, vp0, vs0, vp, vs, dt); 

		
   	 
	 for (iy=0; iy < nyb; iy++) {
	 for (ix=0; ix < nxb; ix++) {
         for (iz=0; iz < nzb; iz++) {
				uu[iy][ix][iz] = upold[iy][ix][iz] + usold[iy][ix][iz];
				vv[iy][ix][iz] = vpold[iy][ix][iz] + vsold[iy][ix][iz];
				ww[iy][ix][iz] = wpold[iy][ix][iz] + wsold[iy][ix][iz];
		}
	}
	}
			
      
   	 uu[isy+nyl][isx+nxl][isz+nbt] += wav[it];
//        ww[isx+nxl][isz+nbt] += wav[it];


		bd3_decay(upold); 
		bd3_decay(upcur); 
		bd3_decay(usold); 
        bd3_decay(uscur); 
		bd3_decay(vpold); 
		bd3_decay(vpcur); 
		bd3_decay(vsold); 
        bd3_decay(vscur); 
		bd3_decay(wpold); 
        bd3_decay(wpcur); 
		bd3_decay(wsold); 
        bd3_decay(wscur); 


        tmp = upold;
        upold = upcur;
        upcur = tmp;

        tmp = vpold;
        vpold = vpcur;
        vpcur = tmp;

        tmp = wpold;
        wpold = wpcur;
        wpcur = tmp;

        tmp = usold;
        usold = uscur;
        uscur = tmp;

        tmp = vsold;
        vsold = vscur;
        vscur = tmp;

        tmp = wsold;
        wsold = wscur;
        wscur = tmp;

		if (it%snap==0) {
		sf_floatwrite(upcur[0][0], nxb*nyb*nzb, fwavup);
		sf_floatwrite(vpcur[0][0], nxb*nyb*nzb, fwavvp);
		sf_floatwrite(wpcur[0][0], nxb*nyb*nzb, fwavwp);
		sf_floatwrite(uscur[0][0], nxb*nyb*nzb, fwavus);
		sf_floatwrite(vscur[0][0], nxb*nyb*nzb, fwavvs);
		sf_floatwrite(wscur[0][0], nxb*nyb*nzb, fwavws);
		sf_floatwrite(uu[0][0], nxb*nyb*nzb, fwavu);
		sf_floatwrite(vv[0][0], nxb*nyb*nzb, fwavv);
		sf_floatwrite(ww[0][0], nxb*nyb*nzb, fwavw);
		}

    }
	
    psstep3e_close();
    bd3_close();

    free(**aa);
    free(*aa);
    free(aa);
    
	free(**vp);     
	free(**vtmp);
    free(**upcur);     
    free(**upold);     
    free(**vpcur);     
    free(**vpold);     
    free(**wpcur);     
    free(**wpold);     
    free(**uscur);     
    free(**usold);     
    free(**vscur);     
    free(**vsold);     
    free(**wscur);     
    free(**wsold);     
    
	free(*vp);     
	free(*vtmp);
    free(*upcur);     
    free(*upold);     
    free(*vpcur);     
    free(*vpold);     
    free(*wpcur);     
    free(*wpold);     
    free(*uscur);     
    free(*usold);     
    free(*vscur);     
    free(*vsold);     
    free(*wscur);     
    free(*wsold);     
    
	free(vp);     
	free(vtmp);
    free(upcur);     
    free(upold);     
    free(vpcur);     
    free(vpold);     
    free(wpcur);     
    free(wpold);     
    free(uscur);     
    free(usold);     
    free(vscur);     
    free(vsold);     
    free(wscur);     
    free(wsold);     
    
fprintf(stderr, "Done\n");
    exit(0); 
}           
           
