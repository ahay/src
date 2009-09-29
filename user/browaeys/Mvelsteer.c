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
    int i,ivx,ivy,it;
    int ixm,iym;

    int nt,nx,ny;
    int nvx,nvy;
    int iypi,iyps,nlive;

    float dt,ot,dx,ox,dy,oy;           
    float dvx,ovx,dvy,ovy;
    float xp,yp,vx,vy;
    float beam,t,x,y;
    float tx,ty;

    float ***data, **semb;

    sf_axis avx,avy,ay;
    sf_file in,out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* input file */
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

    /* output file */
    sf_settype(out,SF_FLOAT);

    if (!sf_getint("nvx",&nvx)) sf_error("Need nvx=");
    /* number of vx values */
    if (!sf_getfloat("dvx",&dvx)) sf_error("Need dvx=");
    /* vx sampling */
    if (!sf_getfloat("ovx",&ovx)) sf_error("Need ovx=");
    /* vx origin */

    if (!sf_getint("nvy",&nvy)) sf_error("Need nvy=");
    /* number of vy values */
    if (!sf_getfloat("dvy",&dvy)) sf_error("Need dvy=");
    /* vy sampling */
    if (!sf_getfloat("ovy",&ovy)) sf_error("Need ovy=");
    /* vy origin */

    if (!sf_getint("iypi",&iypi)) iypi = 0;
    /* first depth position where velocity steering is computed */

    if (!sf_getint("iyps",&iyps)) iyps = ny-1;  
    /* last depth position where velocity steering is computed */
    iyps = SF_MIN(iyps,ny-1);

    avx = sf_maxa(nvx,ovx,dvx);
    sf_oaxa(out,avx,1);

    avy = sf_maxa(nvy,ovy,dvy);
    sf_oaxa(out,avy,2);

    ay = sf_maxa(iyps-iypi+1,oy+iypi*dy,dy);
    sf_oaxa(out,ay,3);

    sf_putstring(out,"label1","vx");
    sf_putstring(out,"label2","vy");
    sf_putstring(out,"label3","y");

    if (!sf_getfloat("xp",&xp)) sf_error("Need xp=");
    /* lateral position where velocity steering is computed */

    /* memory allocations */
    data = sf_floatalloc3(nt,nx,ny);
    semb = sf_floatalloc2(nvx,nvy);

    sf_floatread(data[0][0],nt*nx*ny,in);

    for (i = iypi; i < iyps + 1; i++) {

	yp = oy + i*dy;

	for (ivy = 0; ivy < nvy; ivy++) {
	    for (ivx = 0; ivx < nvx; ivx++) {
		semb[ivy][ivx] = 0.0;
	    }
	}

	for (ivy = 0; ivy < nvy; ivy++) {

	    vy = ovy + ivy*dvy;

	    for (ivx = 0; ivx < nvx; ivx++) {

		vx = ovx + ivx*dvx;

		beam = 0.0;
		nlive = 0;

		for (it = 0; it < nt; it++) {

		    t = it*dt;

                    /* shifted position at current time */
		    x = xp - vx*t;
		    y = yp - vy*t;

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

		} /* it */

                /* normalize stack */
		if (nlive > 0) semb[ivy][ivx] = beam/nlive;

	    } /* vx */

	} /* vy */

	sf_floatwrite(semb[0],nvx*nvy,out);
	sf_warning("i = %d", i);

    } /* y */

    exit(0);
}
