/* 2-D opam for elastic wave modeling and vector field decompostion */
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
#include "opamstep2e.h"
#include "srcsm2.h"
#include "fft2.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nz, nt, ix, iz, it, nbt, nbb, nxl, nxr,  nxb, nzb, isx, isz;
    float dt, dx, dz, o1, o2;
    float **upold,  **upcur, **wpold, **wpcur, *wav;
    float **usold, **uscur, **wsold, **wscur, **tmp;
    float **uu,**ww;
    float  **vp, vp0, **vs, vs0, **vtmp ;
    float ***aa, czt, czb, cxl, cxr; /* top, bottom, left, right */
    float **Cuxxl, **Cuxxr, **Cwxzpl, **Cwxzpr,  **Cuxzpl, **Cuxzpr, **Cwzzl, **Cwzzr, **Cuzzl, **Cuzzr, **Cwxzsl, **Cwxzsr, **Cwxxl, **Cwxxr, **Cuxzsl, **Cuxzsr;
    sf_file Guxxl, Guxxr, Gwxzpl, Gwxzpr,  Guxzpl, Guxzpr,  Gwzzl, Gwzzr, Guzzl, Guzzr, Gwxzsl, Gwxzsr, Gwxxl, Gwxxr, Guxzsl, Guxzsr;
    sf_file fwavup, fwavwp, fwavus, fwavws, fwavu, fwavw, fvelp, fvels, fsource ;
    int opt, snap, nsnap;    /* optimal padding */
    int nkxz,nkxx,nkzz,nxzb,nxzb2,n2,m1, m2, m3, m4, m5, m6, m7, m8;
    int pad1;
    bool cmplx;

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

    Guxxl  = sf_input("Guxxl");
    Guxxr  = sf_input("Guxxr");
    Gwxzpl  = sf_input("Gwxzpl");
    Gwxzpr  = sf_input("Gwxzpr");
    Guxzpl  = sf_input("Guxzpl");
    Guxzpr  = sf_input("Guxzpr");
    Gwzzl  = sf_input("Gwzzl");
    Gwzzr  = sf_input("Gwzzr");
    Guzzl  = sf_input("Guzzl");
    Guzzr  = sf_input("Guzzr");
    Gwxzsl  = sf_input("Gwxzsl");
    Gwxzsr  = sf_input("Gwxzsr");
    Gwxxl  = sf_input("Gwxxl");
    Gwxxr  = sf_input("Gwxxr");
    Guxzsl  = sf_input("Guxzsl");
    Guxzsr  = sf_input("Guxzsr");

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
    if (!sf_getbool("cmplx",&cmplx)) cmplx=false; /* use complex FFT */
    if (!sf_getint("pad1",&pad1)) pad1=1; /* padding factor on the first axis */
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

    nkxz=fft2_init(cmplx,pad1,nzb,nxb,&nkzz,&nkxx);
    nxzb = nxb*nzb;
    nxzb2 = nkxx * nkzz;

    /* check all files */

    if (!sf_histint(Guxxl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Guxxl,"n2",&m1))  sf_error("Need n2= in left");
    
    if (!sf_histint(Guxxr,"n1",&n2) || n2 != m1) sf_error("Need n1=%d in right",m1);
    if (!sf_histint(Guxxr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cuxxl = sf_floatalloc2(nxzb,m1);
    Cuxxr = sf_floatalloc2(m1,nkxz);

    sf_floatread(Cuxxl[0],nxzb*m1,Guxxl);
    sf_floatread(Cuxxr[0],m1*nkxz,Guxxr);

    if (!sf_histint(Gwxzpl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Gwxzpl,"n2",&m2))  sf_error("Need n2= in left");
    
    if (!sf_histint(Gwxzpr,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(Gwxzpr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cwxzpl = sf_floatalloc2(nxzb,m2);
    Cwxzpr = sf_floatalloc2(m2,nkxz);

    sf_floatread(Cwxzpl[0],nxzb*m2,Gwxzpl);
    sf_floatread(Cwxzpr[0],m2*nkxz,Gwxzpr);

    if (!sf_histint(Guxzpl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Guxzpl,"n2",&m3))  sf_error("Need n2= in left");
    
    if (!sf_histint(Guxzpr,"n1",&n2) || n2 != m3) sf_error("Need n1=%d in right",m3);
    if (!sf_histint(Guxzpr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cuxzpl = sf_floatalloc2(nxzb,m3);
    Cuxzpr = sf_floatalloc2(m3,nkxz);

    sf_floatread(Cuxzpl[0],nxzb*m3,Guxzpl);
    sf_floatread(Cuxzpr[0],m3*nkxz,Guxzpr);

    if (!sf_histint(Gwzzl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Gwzzl,"n2",&m4))  sf_error("Need n2= in left");
    
    if (!sf_histint(Gwzzr,"n1",&n2) || n2 != m4) sf_error("Need n1=%d in right",m4);
    if (!sf_histint(Gwzzr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cwzzl = sf_floatalloc2(nxzb,m4);
    Cwzzr = sf_floatalloc2(m4,nkxz);

    sf_floatread(Cwzzl[0],nxzb*m4,Gwzzl);
    sf_floatread(Cwzzr[0],m4*nkxz,Gwzzr);

    if (!sf_histint(Guzzl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Guzzl,"n2",&m5))  sf_error("Need n2= in left");
    
    if (!sf_histint(Guzzr,"n1",&n2) || n2 != m5) sf_error("Need n1=%d in right",m5);
    if (!sf_histint(Guzzr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cuzzl = sf_floatalloc2(nxzb,m5);
    Cuzzr = sf_floatalloc2(m5,nkxz);

    sf_floatread(Cuzzl[0],nxzb*m5,Guzzl);
    sf_floatread(Cuzzr[0],m5*nkxz,Guzzr);

    if (!sf_histint(Gwxzsl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Gwxzsl,"n2",&m6))  sf_error("Need n2= in left");
    
    if (!sf_histint(Gwxzsr,"n1",&n2) || n2 != m6) sf_error("Need n1=%d in right",m6);
    if (!sf_histint(Gwxzsr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cwxzsl = sf_floatalloc2(nxzb,m6);
    Cwxzsr = sf_floatalloc2(m6,nkxz);

    sf_floatread(Cwxzsl[0],nxzb*m6,Gwxzsl);
    sf_floatread(Cwxzsr[0],m6*nkxz,Gwxzsr);

    if (!sf_histint(Gwxxl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Gwxxl,"n2",&m7))  sf_error("Need n2= in left");
    
    if (!sf_histint(Gwxxr,"n1",&n2) || n2 != m7) sf_error("Need n1=%d in right",m7);
    if (!sf_histint(Gwxxr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cwxxl = sf_floatalloc2(nxzb,m7);
    Cwxxr = sf_floatalloc2(m7,nkxz);

    sf_floatread(Cwxxl[0],nxzb*m7,Gwxxl);
    sf_floatread(Cwxxr[0],m7*nkxz,Gwxxr);

    if (!sf_histint(Guxzsl,"n1",&n2) || n2 != nxzb) sf_error("Need n1=%d in left",nxzb);
    if (!sf_histint(Guxzsl,"n2",&m8))  sf_error("Need n2= in left");
    
    if (!sf_histint(Guxzsr,"n1",&n2) || n2 != m8) sf_error("Need n1=%d in right",m8);
    if (!sf_histint(Guxzsr,"n2",&n2) || n2 != nkxz) sf_error("Need n2=%d in right",nkxz);
    Cuxzsl = sf_floatalloc2(nxzb,m8);
    Cuxzsr = sf_floatalloc2(m8,nkxz);

    sf_floatread(Cuxzsl[0],nxzb*m8,Guxzsl);
    sf_floatread(Cuxzsr[0],m8*nkxz,Guxzsr);

//fprintf(stderr, "I'm here...too\n");

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

    lowrankelastic_init2(nzb, nxb, nkxz, nkzz, nkxx, m1, m2, m3, m4, m5, m6, m7, m8,  nxzb2, Cuxxl, Cuxxr, Cwxzpl, Cwxzpr, Cuxzpl, Cuxzpr, Cwzzl, Cwzzr, Cuzzl, Cuzzr, Cwxzsl, Cwxzsr, Cwxxl, Cwxxr, Cuxzsl, Cuxzsr);

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

    srcsm_init(dz, dx);
 
    /* propagation in time */
//	nkxz=fft2_init(true,opt,1,nzb,nxb,&nkzz,&nkxx);
    for (it=0; it < nt; it++) {
	fprintf(stderr, "\b\b\b\b\b%d", it);

	//     uu[isx+nxl][isz+nbt] += wav[it];
	//     ww[isx+nxl][isz+nbt] += wav[it];

        opamstep2e(upold, upcur, wpold, wpcur, usold, uscur, wsold, wscur, uu, ww, nzb, nxb, dz, dx, vp, vs, dt); 
   	 
	for (ix=0; ix < nxb; ix++) {
	    for (iz=0; iz < nzb; iz++) {
		uu[ix][iz] = upold[ix][iz] + usold[ix][iz];
		ww[ix][iz] = wpold[ix][iz] + wsold[ix][iz];
	    }
	}
			
      
   	uu[isx+nxl][isz+nbt] += wav[it];
//      ww[isx+nxl][isz+nbt] += wav[it];
//	 source_smooth(uu, isz+nbt, isx+nxl, wav[it]);
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
    lowrankelastic_close2();
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
           
