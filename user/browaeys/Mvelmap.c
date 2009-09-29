/* 2-D mapping from moving-object velocity to plane-wave slowness */
/*
  Copyright (C) 2009 University of Texas at Austin

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

static const float eps = 1.e-5;

int main(int argc, char* argv[])
{
    int isx,isy,it;                /* loop counters */
    int ivxm,ivym;                 /* cell indices for bilinear interpolation */

    int nsx,nsy;                   /* number of slowness in x and y */
    int nvx,nvy;                   /* number of velocity in x and y */
    int nz,i;

    int nt;                        /* number of line parameters */

    float dt,ot;                   /* increment of parameter, starting parameter */
    float dvx,ovx;                 /* increment of velocity in x, starting velocity */
    float dvy,ovy;                 /* increment of velocity in y, starting velocity */
    float oz,dz;

    float vxmax,vymax;
    float osx,osy,dsx,dsy;
    float sx,sy,vx,vy;
    float t,c,s,dn,line,tvx,tvy;

    float **vv;                    /* value on velocity grid vx-vy */
    float **sv;                    /* value on slowness grid sx-sy */

    sf_axis asx,asy;
    sf_file in,out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nvx)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dvx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&ovx)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&nvy)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dvy)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&ovy)) sf_error("No o2= in input");

    if (!sf_histint(in,"n3",&nz)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dz)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&oz)) oz=0.;

    if ((nvx != nvy) || (dvx != dvy) || (ovx != ovy)) sf_error("Need centered squared grid for vx-vy plane");

    vxmax = ovx + (nvx-1)*dvx;
    vymax = ovy + (nvy-1)*dvy;

    /* slowness grid */
    dsx = 1./vxmax;
    dsy = 1./vymax;

    if (!sf_getfloat("osx",&osx)) osx = -0.5/dvx;
    if (!sf_getfloat("osy",&osy)) osy = -0.5/dvy;

    if ((osx != osy) || (osx >= 0.0) || (osy >= 0.0)) sf_error("Need centered squared grid for sx-sy plane");

    nsx = 2*floorf(fabs(osx)/dsx);
    nsy = 2*floorf(fabs(osy)/dsy);

    /* slowness grid do not contain zero */
    osx = osx + 0.5*dsx;
    osy = osy + 0.5*dsy;

    /* line trigonometric parametrization */
    if (!sf_getint("nt",&nt)) nt=360;
    /* number of line parameter for integration in [0,180]. */
    if (!sf_getfloat("dt",&dt)) dt=0.5;
    /* line parameter increment */
    if (!sf_getfloat("ot",&ot)) ot=0.;
    /* line parameter origin */

    /* output file parameters */ 
    asx = sf_maxa(nsx,osx,dsx);
    sf_oaxa(out,asx,1);

    asy = sf_maxa(nsy,osy,dsy);
    sf_oaxa(out,asy,2);
    
    sf_putstring(out,"label1","sx");
    sf_putstring(out,"label2","sy");

     /* memory allocations */
    vv = sf_floatalloc2(nvx,nvy);
    sv = sf_floatalloc2(nsx,nsy);

    for (i = 0; i < nz; i++) {

        /* clear array and read data */
	for (isy = 0; isy < nsy; isy++) {
	    for (isx = 0; isx < nsx; isx++) {
                sv[isy][isx] = 0.;
	    }
	}

	sf_floatread(vv[0],nvx*nvy,in);

	for (isy = 0; isy < nsy; isy++) {

	    sy = osy + isy*dsy;

	    for (isx = 0; isx < nsx; isx++) {

		sx = osx + isx*dsx;

		line = 0.0;

                /* loop on velocity in the line */
		for (it = 0; it < nt; it++) {

                    /* line parameter */
		    t = (ot + it*dt)*SF_PI/180.0;

		    c = cos(t);
		    s = sin(t);

		    dn = 1.0/(c*sx + s*sy + eps);
		    vx = c*dn;
		    vy = s*dn;

                    /* bilinear interpolation */
		    ivxm = floorf((vx-ovx)/dvx);
		    ivym = floorf((vy-ovy)/dvy);

		    if ( (ivxm >= 0) && (ivym >= 0) && (ivxm < (nvx-1)) && (ivym < (nvy-1)) ) {

			tvx = (vx-ivxm*dvx-ovx)/dvx;
			tvy = (vy-ivym*dvy-ovy)/dvy;

                        /* add contribution from this velocity into line sum */
			line += vv[ivym][ivxm]*(1.-tvx)*(1.-tvy);
			line += vv[ivym][ivxm+1]*tvx*(1.-tvy);
			line += vv[ivym+1][ivxm+1]*tvx*tvy;
			line += vv[ivym+1][ivxm]*(1.-tvx)*tvy;
		    }

		}

		sv[isy][isx] = line/nt;

	    }
	}

        /* output on slowness grid */
	sf_floatwrite(sv[0],nsx*nsy,out);
	sf_warning("i = %d", i);

    }

    exit(0);
}
