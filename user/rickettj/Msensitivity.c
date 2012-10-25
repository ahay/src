/* Traveltime sensitivity kernels. */
/*
  Copyright (C) 2000 The Board of Trustees of Stanford University
  
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

#include "scattering.h"
#include "raytrace.h"
#include "green.h"

#define ONLY_ONE 0
#define EPS 1.0e-6

static void MinimumTravelTime(float xx,float xy,float xz,float *output);
static void TravelTime(float xx,float xy,float xz,float *output);

static float sx,sy,sz,rx,ry,rz;

int main (int argc, char* argv[]) 
{
    int ix,iy,iz;
    int mode;

    float xx,xy,xz;

    float *output;
    double tt1,ampl;
    int ipx,ipy,ipz;
    int nx,ny,nz,/* nv, */ nt,nw;
    float ox,oy,oz,dx,dy,dz,dv;
    float v0,vgrad;
    float t1,t2,dt,ot,dw,ow;
    float px,py,pz;
    sf_file out;

    sf_init(argc, argv);

    if (!sf_getint("nx",&nx)) nx=100;
    if (!sf_getint("ny",&ny)) ny=1;
    if (!sf_getint("nz",&nz)) nz=100; /* dimensions */
    if (!sf_getfloat("ox",&ox)) ox=0.;
    if (!sf_getfloat("oy",&oy)) oy=0.;
    if (!sf_getfloat("oz",&oz)) oz=0.; /* grid origin */
    if (!sf_getfloat("dx",&dx)) dx=0.02;
    if (!sf_getfloat("dy",&dy)) dy=0.02;
    if (!sf_getfloat("dz",&dz)) dz=0.02; /* grid spacing */
    dv=dx*dy*dz;
/*    nv=nx*ny*nz; */

    if (!sf_getfloat("v0",   &v0   )) v0=1.5;   /* surface velocity */
    if (!sf_getfloat("vgrad",&vgrad)) vgrad=.8; /* velocity gradient */
    GreenInit(v0,vgrad);

    if (!sf_getint("mode",&mode)) mode=6;

    /* source coordinates */
    if (!sf_getfloat("sx",&sx)) sx=0.15;
    if (!sf_getfloat("sy",&sy)) sy=0.;
    if (!sf_getfloat("sz",&sz)) sz=0.15; /* source */

    /* receiver coordinates */
    if (!sf_getfloat("rx",&rx)) rx=1.6;
    if (!sf_getfloat("ry",&ry)) ry=0.;
    if (!sf_getfloat("rz",&rz)) rz=1.605; /* receiver */

    /* other parameters */
    if (!sf_getfloat("t1",&t1)) t1=0.8;
    if (!sf_getfloat("t2",&t2)) t2=1.2;
    if (!sf_getfloat("dt",&dt)) dt=0.008;
    ot=t1; nt=1+(t2-t1)/dt;

    nw=nt;
    dw=2*SF_PI/(2*nt*dt);
    ow=0.;

    ipx=1;  ipy=2;  ipz=24;
    px=ipx*dx+ox;
    py=ipy*dy+oy; 
    pz=ipz*dz+oz;
    px=0.7;py=0.;pz=0.7+dz/2;

    if ((mode==2)||(mode==5)||(mode==6)) 
	ScatteringInit(dv,
		       nt,ot,dt,
		       nw,ow,dw,
		       sx,sy,sz,
		       rx,ry,rz,
		       px,py,pz);

    out = sf_output("out");
    sf_setformat(out,"native_float");
    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);
    sf_putint(out,"n3",nz);
    sf_putfloat(out,"o1",ox);
    sf_putfloat(out,"o2",oy);
    sf_putfloat(out,"o3",oz);
    sf_putfloat(out,"d1",dx);
    sf_putfloat(out,"d2",dy);
    sf_putfloat(out,"d3",dz);

    output=sf_floatalloc(nx);

    GreenTtAmp(rx,ry,rz,sx,sy,sz,&tt1,&ampl);

    sf_warning("Background traveltime: %f",tt1);

    if (mode==3) { 
	RayTrace(nx,ny,nz,
		 ox,oy,oz,
		 dx,dy,dz,
		 rx,ry,rz,
		 sx,sy,sz,
		 v0,vgrad,
		 out); 
	exit(0); 
    }

    for     (iz=0; iz<nz; iz++) { 
	for   (iy=0; iy<ny; iy++) {  
	    for (ix=0; ix<nx; ix++) {  
		xz=iz*dz+oz;	xy=iy*dy+oy;	xx=ix*dx+ox;


		if (((ipx==ix)&&(ipy==iy)&&(ipz==iz))||(ONLY_ONE!=1)) {

		    switch(mode) {
			case 1:
			    MinimumTravelTime(xx,xy,xz,(output+ix)); 
			    break;
			case 2: 
			    BornSensitivity(xx,xy,xz,(output+ix));
			    break;
			case 6:   
			    RytovSensitivity(xx,xy,xz,(output+ix));
			    break;
			case 4:
			    TravelTime(xx,xy,xz,(output+ix));
			    break;
			case 5:
			    BornScatteredField(xx,xy,xz,(output+ix));
			    break;
			default:
			    sf_error("mode %d not implemented",mode);
		    }
		} else {
		    output[ix]=0.;
		}	  
	    }
	    sf_floatwrite(output,nx,out);
	}
    }

    exit(0);
}

static void MinimumTravelTime(float xx,float xy,float xz,float *output)
{
    double tt0,tt1,tt2,ampl;
    GreenTtAmp(rx,ry,rz,sx,sy,sz,&tt0,&ampl);

    GreenTtAmp(rx,ry,rz,xx,xy,xz,&tt1,&ampl);
    GreenTtAmp(xx,xy,xz,sx,sy,sz,&tt2,&ampl);
    tt1+=tt2;

    if (tt1<(tt0+100*EPS)) *output=1.;
    else 	              *output=0.;
    return;
}

static void TravelTime(float xx,float xy,float xz,float *output)
{
    double tt2,ampl;
    GreenTtAmp(xx,xy,xz,sx,sy,sz,&tt2,&ampl);
    *output=(float) tt2;
    return;
}




