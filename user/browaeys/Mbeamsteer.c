/* Beam steering for 2D surface array. */
/*
  Copyright (C) 2008 University of Texas at Austin
  Adapted from Steve Cole, Stanford University, 1995

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

#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)

int main(int argc, char* argv[])
{
    int nt,it;                     /* number of time samples, time counter */
    int nx,ix;                     /* number of receivers in x, receivers counter */
    int ny,iy;                     /* number of receivers in y, receivers counter */
    int nlive;                     /* number of non zero traces */
    int npx,npy;                   /* number of slopes in x, number of slopes in y */

    float dt,ot;                   /* time increment, starting time */
    float dx,ox;                   /* increment in x, starting position */
    float dy,oy;                   /* increment in y, starting position */
    float xref,yref;               /* coordinates where beams are computed */
    float dpx,opx;                 /* slope increment, starting slope in x */
    float dpy,opy;                 /* slope increment, starting slope in y */
    float sum;                     /* test sum for zero traces */

    bool mode;                     /* 2-D slope vectors, or slowness and azimuth */ 
    bool **live;                   /* non zero traces flag */
    float ***data;                 /* 2-D surface seismic data */
    float **semb;                  /* focused slant stack */

    float x,y;                     /* receiver positions */
    float px,py;                   /* slopes in x and y */
    float tshift;                  /* time shift */
    float beam;                    /* beam stack */

    int istart,istop;              /* time integration start and end indices */
    int ipx,ipy;                   /* slopes counters */
    int itshift;                   /* time shift in samples */

    sf_file in,out;
   
    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("mode",&mode)) mode=true;
    /* if n, beams computed as a function of apparent slowness and azimuth angle. */
      
    if (!sf_getfloat("xref",&xref)) sf_error("Need xref=");
    /* x coordinate where beams are computed */

    if (!sf_getfloat("yref",&yref)) sf_error("Need yref=");
    /* y coordinate where beams are computed */

    if (!sf_getint("npx",&npx)) sf_error("Need npx=");
    /* number of px values (if mode=y). */
    if (!sf_getfloat("dpx",&dpx)) sf_error("Need dpx=");
    /* px sampling (if mode=y). */
    if (!sf_getfloat("opx",&opx)) sf_error("Need opx=");
    /* px origin (if mode=y) */

    if (!sf_getint("npy",&npy)) sf_error("Need npy=");
    /* number of py values (if mode=y). */
    if (!sf_getfloat("dpy",&dpy)) sf_error("Need dpy=");
    /* py sampling (if mode=y). */
    if (!sf_getfloat("opy",&opy)) sf_error("Need opy=");
    /* py origin (if mode=y) */  

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&ot)) ot=0.;

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&ox)) ox=0.;

    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&oy)) oy=0.;

    /* output file parameters */ 
    sf_putint(out,"n1",npx);
    sf_putfloat(out,"d1",dpx);
    sf_putfloat(out,"o1",opx);

    sf_putint(out,"n2",npy);
    sf_putfloat(out,"d2",dpy);
    sf_putfloat(out,"o2",opy);

    /* memory allocations */
    data = sf_floatalloc3(nt,nx,ny);
    semb = sf_floatalloc2(npx,npy);
    live = sf_boolalloc2(nx,ny);

    /* clear data array and read the data */
    for (iy = 0; iy < ny; iy++) {
	for (ix = 0; ix < nx; ix++) {
	    for (it = 0; it < nt; it++) {
		data[iy][ix][it] = 0.;
	    }
	}
    }
    sf_floatread(data[0][0],ny*nx*nt,in);

    /* determine if this is a dead or live trace */
    nlive = 0;
    for (iy = 0; iy < ny; iy++) {
	for (ix = 0; ix < nx; ix++) {
	    sum = 0.;
	    for (it = 0; it < nt; it++) {
		sum += data[iy][ix][it];
	    }
	    if (sum == 0.0) {
		live[iy][ix] = false;
	    } else {
		live[iy][ix] = true;
		nlive += 1;
	    }
	}
    }

    /* beam steering */

    /* loop over 2-D slopes */
    for (ipy = 0; ipy < npy; ipy++) {
	py = opy + ipy*dpy;

	for (ipx = 0; ipx < npx; ipx++) {
	    px = opx + ipx*dpx;

            /* clear arrays */
	    semb[ipy][ipx] = 0.;
	    beam = 0.;
	    
	    /* loop over receivers */
	    for (iy = 0; iy < ny; iy++) {

		for (ix = 0; ix < nx; ix++) {

		    if (live[iy][ix]) {

                        /* position compared to reference trace */
			y = iy*dy + oy - yref;
			x = ix*dx + ox - xref;

                        /* compute the necessary time shift to align */
			/* the current trace with the reference trace */
			if (mode) {
			    tshift = px*x + py*y;
			} else {
                            /* py is apparent slowness, px is azimuth */
			    tshift = py*( cosf(SF_PI*px/180.)*x + sinf(SF_PI*px/180.)*y );
			}

                        /* Nearest integer */
			itshift = floorf(tshift/dt);
			if ( (2.*tshift) > ((2*itshift+1)*dt) ) itshift += 1;

                        /* loop over samples of beam */
			istart = MAX(itshift,0);
			istop = MIN(nt+itshift,nt);
			for (it = istart; it < istop; it++){
                            /* sum the shifted input into the beam */
			    beam += data[iy][ix][it-itshift];
			}

		    }
		}
	    }

            /* normalize stack */
	    semb[ipy][ipx] = beam/nlive;

	}
    }

    /* output beam stack */
    sf_floatwrite(semb[0],npy*npx,out);

    exit(0);
}
