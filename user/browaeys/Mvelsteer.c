/* Velocity steering for 2D receivers array. */
/*
  Copyright (C) 2008 University of Texas at Austin

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
#include <float.h>

#include <rsf.h>


int main(int argc, char* argv[])
{
    int nt,it;                     /* number of time samples, time counter */
    int nx,ix;                     /* number of receivers in x, receivers counter */
    int ny,iy;                     /* number of receivers in y, receivers counter */
    int nvx,nvy;                   /* number of velocity in x, number of velocity in y */
    int nref;                      /* number of points where steering is preformed */
    int nlive;                     /* number of summed samples in the beam */

    float dt,ot;                   /* time increment, starting time */
    float dx,ox;                   /* increment in x, starting position */
    float dy,oy;                   /* increment in y, starting position */
    float xref,yref;               /* coordinates where beams are computed */
    float dvx,ovx;                 /* velocity increment, starting velocity in x */
    float dvy,ovy;                 /* velocity increment, starting velocity in y */

    bool mode;                     /* 2-D velocity vectors, or velocity and angle */ 
    float ***data;                 /* 2-D surface seismic data */
    float ***semb;                 /* focused beam stack */

    float t;                       /* time position */
    float vx,vy;                   /* velocity along x and y */
    float x,y;                     /* x,y positions */
    float tx,ty;                   /* x,y position inside cell for interpolation */
    float beam;                    /* beam trajectory stack */

    int ivx,ivy;                   /* velocity counters */
    int ixm,iym;                   /* shifted space positions in samples */
    
    sf_axis avx,avy,aref;
    sf_file in,out;

    sf_init(argc,argv);

    /* Axis label 1:t, 2:x, 3:y */
    in = sf_input("in");

    /* Axis label 1:vx, 2:vy 3:ref*/
    out = sf_output("out");

    if (!sf_getbool("mode",&mode)) mode=true;
    /* if n, beams computed as a function of apparent velocity and angle. */
      
    if (!sf_getfloat("xref",&xref)) sf_error("Need xref=");
    /* x coordinate where beams are computed */

    if (!sf_getfloat("yref",&yref)) sf_error("Need yref=");
    /* y coordinate where beams are computed */

    if (!sf_getint("nvx",&nvx)) sf_error("Need nvx=");
    /* number of vx values (if mode=y). */
    if (!sf_getfloat("dvx",&dvx)) sf_error("Need dvx=");
    /* vx sampling (if mode=y). */
    if (!sf_getfloat("ovx",&ovx)) sf_error("Need ovx=");
    /* vx origin (if mode=y) */

    if (!sf_getint("nvy",&nvy)) sf_error("Need nvy=");
    /* number of vy values (if mode=y). */
    if (!sf_getfloat("dvy",&dvy)) sf_error("Need dvy=");
    /* vy sampling (if mode=y). */
    if (!sf_getfloat("ovy",&ovy)) sf_error("Need ovy=");
    /* vy origin (if mode=y) */  

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&ot)) ot=0.;

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&ox)) ox=0.;

    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&oy)) oy=0.;

    /* number of points where beams are computed. */
    nref =1;

    /* output file parameters */ 
    avx = sf_maxa(nvx,ovx,dvx);
    sf_oaxa(out,avx,1);

    avy = sf_maxa(nvy,ovy,dvy);
    sf_oaxa(out,avy,2);

    aref = sf_maxa(nref,0,1);
    sf_oaxa(out,aref,3);

    sf_putstring(out,"label1","vx");
    sf_putstring(out,"label2","vy");
    sf_putstring(out,"label3","ref");

    /* memory allocations */
    data = sf_floatalloc3(nt,nx,ny);
    semb = sf_floatalloc3(nvx,nvy,nref);

    /* clear data array and read the data */
    for (iy = 0; iy < ny; iy++) {
	for (ix = 0; ix < nx; ix++) {
	    for (it = 0; it < nt; it++) {
		data[iy][ix][it] = 0.;
	    }
	}
    }
    sf_floatread(data[0][0],nt*nx*ny,in);

    /* velocity beam steering */

    /* loop over 2-D velocity */
    for (ivy = 0; ivy < nvy; ivy++) {
	vy = ovy + ivy*dvy;

	for (ivx = 0; ivx < nvx; ivx++) {
	    vx = ovx + ivx*dvx;

            /* clear arrays */
	    semb[0][ivy][ivx] = 0.;
	    beam = 0.;
	    
	    /* loop over time samples */
	    nlive = 0;
	    for (it = 0; it < nt; it++) {

                /* time position */
		t = it*dt;

                /* compute shifted position at current time */
		if (mode) {
		    x = xref - vx*t;
		    y = yref - vy*t;
		} else {
                    /* vy is velocity, vx is dip angle */
		    x = xref - vy*cosf(SF_PI*vx/180.)*t;
		    y = yref - vy*sinf(SF_PI*vx/180.)*t;
		}

                /* bilinear interpolation */
		ixm = floorf((x-ox)/dx);
		iym = floorf((y-oy)/dy);

		if ( (ixm >= 0) && (iym >= 0) && (ixm < (nx-1)) && (iym < (ny-1)) ) {

                       tx = (x-ixm*dx-ox)/dx;
                       ty = (y-iym*dy-oy)/dy;

                       /* sum the directional input at t into the beam */
                       beam += data[iym][ixm][it]*(1.-tx)*(1.-ty);
                       beam += data[iym][ixm+1][it]*tx*(1.-ty);
                       beam += data[iym+1][ixm+1][it]*tx*ty;
                       beam += data[iym+1][ixm][it]*(1.-tx)*ty;
		       nlive += 1;
		}
	    }

            /* normalize stack */
	    semb[0][ivy][ivx] = beam/nlive;

	}
    }

    /* output velocity beam stack */
    sf_floatwrite(semb[0][0],nvx*nvy*1,out);

    exit(0);
}
