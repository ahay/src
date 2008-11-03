/* Convert F(vx,vy) to F(d) using equation vy*cos(d) + vx*sin(d) = 1/s0. */
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

#include <rsf.h>

int main(int argc, char* argv[])
{
    int na;                        /* number of output angles */
    int nt;                        /* number of radial lines or polar angles */
    int npos;                      /* number of frames */
    int nvx,nvy;                   /* number of velocity in x, number of velocity in y */

    float da,oa;                   /* increment in output angle, starting position */
    float dt,ot;                   /* increment in polar angle, starting position */
    float dvx,ovx;                 /* increment of velocity in x, starting velocity */
    float dvy,ovy;                 /* increment of velocity in y, starting velocity */
    float s0;                      /* reference slowness */

    float **v;                     /* data on squared grid vx-vy */
    float **a;                     /* data on output angle grid */

    int ia;                        /* output angle counter */
    int it;                        /* polar angle counter */
    int ivxm,ivym;                 /* bilinear interpolation cell indices */

    float d;                       /* output angle */
    float t;                       /* polar angle defining radial line */
    float vi;                      /* inverse velocity radius */
    float dat;                     /* interpolated data value at vx,vy */
    float vx,vy;                   /* velocity vx and vy */
    float tx,ty;                   /* bilinear inside cell position */

    sf_axis aa,apos;
    sf_file in,out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
      
    if (!sf_getfloat("s0",&s0)) sf_error("Need s0=");
    /* reference slowness */

    if (!sf_getint("na",&na)) sf_error("Need na=");
    /* number of output angle values. */
    if (!sf_getfloat("da",&da)) sf_error("Need da=");
    /* output angle sampling. */
    if (!sf_getfloat("oa",&oa)) sf_error("Need oa=");
    /* output angle origin */

    if (!sf_getint("nt",&nt)) nt=360;
    /* number of polar angles for integration. */
    if (!sf_getfloat("dt",&dt)) dt=1.;
    /* polar angle sampling. */
    if (!sf_getfloat("ot",&ot)) ot=0.;
    /* polar angle origin */

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nvx)) sf_error("No nvx= in input");
    if (!sf_histfloat(in,"d1",&dvx)) sf_error("No dvx= in input");
    if (!sf_histfloat(in,"o1",&ovx)) sf_error("No ovx= in input");

    if (!sf_histint(in,"n2",&nvy)) sf_error("No nvy= in input");
    if (!sf_histfloat(in,"d2",&dvy)) sf_error("No dvy= in input");
    if (!sf_histfloat(in,"o2",&ovy)) sf_error("No ovy= in input");

    if ((nvx != nvy) || (dvx != dvy)) sf_error("Need squared grid for vx-vy plane");

    /* number of vx-vy frames for computation. */
    npos = 1;

    /* output file parameters */ 
    aa = sf_maxa(na,oa,da);
    sf_oaxa(out,aa,1);

    apos = sf_maxa(npos,0,1);
    sf_oaxa(out,apos,2);
    
    sf_putstring(out,"label1","angle");
    sf_putstring(out,"label2","pos");

    /* memory allocations */
    v = sf_floatalloc2(nvx,nvy);
    a = sf_floatalloc2(na,npos);

    /* read data in velocity array */
    sf_floatread(v[0],nvx*nvy,in);

    /* initialize */
    for (ia = 0; ia < na; ia++) {
	a[0][ia] = 0.;
    }

    /* loop on polar angle directions */
    for (it = 0; it < nt; it++) {
	t = (ot + it*dt)/180.*SF_PI;

        /* loop on output angle */
	for (ia = 0; ia < na; ia++) {
	    d = (oa + ia*da)/180.*SF_PI;

            /* inverse velocity radius equals scaled cos(d-t) on radial line */
	    vi = s0*cosf(d-t);

	    if (vi != 0.0) {

		vx = sinf(t)/vi;
		vy = cosf(t)/vi;
 
                /* bilinear interpolation */
		ivxm = floorf((vx-ovx)/dvx);
		ivym = floorf((vy-ovy)/dvy);

		if ( (ivxm >= 0) && (ivym >= 0) && (ivxm < (nvx-1)) && (ivym < (nvy-1)) ) {

		    tx = (vx-ivxm*dvx-ovx)/dvx;
		    ty = (vy-ivym*dvy-ovy)/dvy;

                    /* calculate interpolated value at vx,vy */
		    dat = 0.0;
		    dat += v[ivym][ivxm]*(1.-tx)*(1.-ty);
		    dat += v[ivym][ivxm+1]*tx*(1.-ty);
		    dat += v[ivym+1][ivxm+1]*tx*ty;
		    dat += v[ivym+1][ivxm]*(1.-tx)*ty;

                    /* sum */
		    a[0][ia] += dat;
		}
	    }
	}
    }
	
    /* output on angle grid */
    sf_floatwrite(a[0],na*npos,out);

    exit(0);
}
