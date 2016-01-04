/* Forward and inverse normal moveout with interval velocity. */
/*
  Copyright (C) 2004 University of Texas at Austin

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

#include <math.h>
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "orthovelo.h"
#include "warp3.h"

int main (int argc, char* argv[])
{
    bool inv,cij;
    int it,ipx, ipy, nt,nx, ny, npx, npy, ileft, nleft;
    float dt, t0, px, px0, px2, py, py0, py2 , v, vx, vy, vz, vgp, sx, st, sy, ft, dpx, dpy, eps, x0, dx, y0, dy;
    float c11, c22, c33, c44, c55, c66, c12, c13, c23;
    float ***tp, *vel, *c_11, *c_22, *c_33, *c_44, *c_55, *c_66, *c_12, *c_13, *c_23, ***txy, ***t, ***x, ***y;
    sf_file inp, out;
    sf_file C11,C22,C33,C44,C55,C66,C12,C13,C23, velocity;

    sf_init (argc,argv);

    inp = sf_input("in");
    out = sf_output("out");


    if (!sf_getbool("cij",&cij)) cij=false;
    if (cij) { /* Orthorhombic */
        C11 = sf_input("c11");
        C22 = sf_input("c22");
        C33 = sf_input("c33");
        C44 = sf_input("c44");
        C55 = sf_input("c55");
        C66 = sf_input("c66");
        C12 = sf_input("c12");
        C13 = sf_input("c13");
        C23 = sf_input("c23");
    } else { /* Isotropic */
        velocity = sf_input("velocity");
   }


    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_getbool("inv",&inv)) inv=false;

    if (inv) { /* tp -> tx */
    	if (!sf_histint(inp,"n2",&npx)) sf_error("No n2= in input");
    	if (!sf_histfloat(inp,"d2",&dpx)) sf_error("No d2= in input");
    	if (!sf_histfloat(inp,"o2",&px0)) sf_error("No o2= in input");
    	
    	if (!sf_histint(inp,"n3",&npy)) sf_error("No n3= in input");
    	if (!sf_histfloat(inp,"d3",&dpy)) sf_error("No d3= in input");
    	if (!sf_histfloat(inp,"o3",&py0)) sf_error("No o3= in input");
    	
    	if (!sf_getint("nx",&nx)) sf_error("Need nx=");
    	/* x offset samples */
    	if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
    	/* x offset sampling */
    	if (!sf_getfloat("x0",&x0)) x0=0.;
    	/* x first offset */
    	
    	if (!sf_getint("ny",&ny)) sf_error("Need ny=");
    	/* y offset samples */
    	if (!sf_getfloat("dy",&dy)) sf_error("Need dy=");
    	/* y offset sampling */
    	if (!sf_getfloat("y0",&y0)) y0=0.;
    	/* y first offset */
    	
    	sf_putint(out,"n2",nx);
    	sf_putfloat(out,"d2",dx);
    	sf_putfloat(out,"o2",x0);
    	sf_putint(out,"n3",ny);
    	sf_putfloat(out,"d3",dy);
    	sf_putfloat(out,"o3",y0);
    	
    } else {  /* tx -> tp */
    	if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
    	if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in input");
    	if (!sf_histfloat(inp,"o2",&x0)) sf_error("No o2= in input");
    	
    	if (!sf_histint(inp,"n3",&ny)) sf_error("No n3= in input");
    	if (!sf_histfloat(inp,"d3",&dy)) sf_error("No d3= in input");
    	if (!sf_histfloat(inp,"o3",&y0)) sf_error("No o3 in input");
    	
    	if (!sf_getint("npx",&npx)) sf_error("Need npx=");
    	/* x slope samples */
    	if (!sf_getfloat("dpx",&dpx)) sf_error("Need dpx=");
    	/* x slope sampling */
    	if (!sf_getfloat("px0",&px0)) px0=0.;
    	/* x first slope */
    	
    	if (!sf_getint("npy",&npy)) sf_error("Need npy=");
    	/* y slope samples */
    	if (!sf_getfloat("dpy",&dpy)) sf_error("Need dpy=");
    	/* y slope sampling */
    	if (!sf_getfloat("py0",&py0)) py0=0.;
    	/* y first slope */
    	
    	sf_putint(out,"n2",npx);
    	sf_putfloat(out,"d2",dpx);
    	sf_putfloat(out,"o2",px0);
    	
    	sf_putint(out,"n3",npy);
    	sf_putfloat(out,"d3",dpy);
    	sf_putfloat(out,"o3",py0);
    }

    nleft = sf_leftsize(inp,3);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    warp3_init(nt, t0, dt,
	       nx, x0, dx,
	       ny, y0, dy /* output grid */,
	       nt, npx, npy     /* input  grid */,
	       eps     /* regularization */);

    tp = sf_floatalloc3(nt,npx,npy);
    t = sf_floatalloc3(nt,npx,npy);
    x = sf_floatalloc3(nt,npx,npy);
    y = sf_floatalloc3(nt,npx,npy);

    txy   = sf_floatalloc3(nt,nx,ny);
    
    if (cij) { /* Orthorhombic */
        c_11 = sf_floatalloc(nt);
        c_22 = sf_floatalloc(nt);
        c_33 = sf_floatalloc(nt);
        c_44 = sf_floatalloc(nt);
        c_55 = sf_floatalloc(nt);
        c_66 = sf_floatalloc(nt);
        c_12 = sf_floatalloc(nt);
        c_13 = sf_floatalloc(nt);
        c_23 = sf_floatalloc(nt);
    } else { /* Isotropic */
        vel = sf_floatalloc(nt);
   }


    for (ileft = 0; ileft < nleft; ileft++) {

	    if (cij) { /* Orthorhombic */
            sf_floatread (c_11,nt,C11);
            sf_floatread (c_22,nt,C22);
            sf_floatread (c_33,nt,C33);
            sf_floatread (c_44,nt,C44);
            sf_floatread (c_55,nt,C55);
            sf_floatread (c_66,nt,C66);
            sf_floatread (c_12,nt,C12);
            sf_floatread (c_13,nt,C13);
            sf_floatread (c_23,nt,C23);
        } else { /* Isotropic */
	        sf_floatread (vel,nt,velocity);
       }
       
/*  #pragma omp parallel for default(shared)				\
    private(ipy,py,ipx,px,st,sx,sy,it,c11,c22,c33,c44,c55,c66,c12,c13,c23,vx,vy,vz,vgp,v,py2,px2,ft) */
    for (ipy = 0; ipy < npy; ipy++) {
	    py = py0 + ipy*dpy;
	    py2 = py*py;
	    
    	for (ipx = 0; ipx < npx; ipx++) {
    	    px = px0 + ipx*dpx;
    	    px2= px*px;

    	    st = 0.;
    	    sx = 0.;
    	    sy = 0.;
    	    for (it=0; it < nt; it++) {
	    	    if (cij) { /* Orthorhombic */
                    c11 = c_11[it]; c22 = c_22[it]; c33 = c_33[it]; c44 = c_44[it]; c55 = c_55[it]; 
                    c66 = c_66[it]; c12 = c_12[it]; c13 = c_13[it]; c23 = c_23[it];  
                    vx = xPvelo(c11,c22,c33,c44,c55,c66,c12,c13,c23,px,py); /* group velocity computation*/
                    vy = yPvelo(c11,c22,c33,c44,c55,c66,c12,c13,c23,px,py);
                    vz = zPvelo(c11,c22,c33,c44,c55,c66,c12,c13,c23,px,py);
                    vgp = sqrtf(vx*vx+vy*vy+vz*vz);
            		if (it==0) {
            		    st = t0*(vgp/vz)/dt;
            		    sx = vx*st;
            		    sy = vy*st;
            		}

            		t[ipy][ipx][it] = st*dt;
            		x[ipy][ipx][it] = sx*dt;
            		y[ipy][ipx][it] = sy*dt;
            		
            		st += (vgp/vz);
            		sx += vx*(vgp/vz);
            		sy += vy*(vgp/vz);
            		
                } else { /* Isotropic */
                
            		v = vel[it];
            		v *= v;

            		ft = 1.0-(px2+py2)*v;

            		if (ft < 0.) { /* Ignore from consideration */
            		    for (it=0; it < nt; it++) {
                			t[ipy][ipx][it]=t0-(nt+10)*dt;
                			x[ipy][ipx][it]=x0-(nx+10)*dx;
                			y[ipy][ipx][it]=y0-(ny+10)*dy;
            		    }
            		    break;
            		}

            		ft = sqrtf(1.0/ft);

            		if (it==0) {
            		    st = t0*ft/dt;
            		    sx = px*v*st;
            		    sy = py*v*st;
            		}

            		t[ipy][ipx][it] = st*dt;
            		x[ipy][ipx][it] = sx*dt;
            		y[ipy][ipx][it] = sy*dt;
            		
            		st += ft;
            		sx += px*v*ft;
            		sy += py*v*ft;
               }
    	    }
    	}
	}
	


	if (inv) {
	    sf_floatread (tp[0][0],nt*npx*npy,inp);	    
	    warp3(tp,t,x,y,txy);	    
	    sf_floatwrite (txy[0][0],nt*nx*ny,out);
	} else {
	    sf_floatread (txy[0][0],nt*nx*ny,inp);
	    fwarp3(txy,t,x,y,tp);
	    sf_floatwrite (tp[0][0],nt*npx*npy,out);
	}
    }

/*	 for (ipy = 0; ipy < npy; ipy++) {*/
/*    	for (ipx = 0; ipx < npx; ipx++) {*/
/*    	    for (it=0; it < nt; it++) {*/
/*            		x[ipy][ipx][it] = abs(x[ipy][ipx][it])-abs(y[ipy][ipx][it]);*/
/*               }*/
/*    	    }*/
/*    	}*/
sf_file testtime, testdiff1,testdiff2;
testtime = sf_output("time");
testdiff1 = sf_output("diff1");
testdiff2 = sf_output("diff2");
sf_floatwrite(t[0][0],nt*npx*npy,testtime);
sf_floatwrite(x[0][0],nt*npx*npy,testdiff1);
sf_floatwrite(y[0][0],nt*npx*npy,testdiff2);



    exit (0);
}
