/* 2-D slowness vector to angle transformation. */
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

int main (int argc, char* argv[])
{
    int k,ia,i;

    int nsx,nsy,nz;
    int na,nr,nlive,isxm,isym;

    float osx,dsx,osy,dsy,oz,dz;
    float oa,da,orig,dr;
    float a,sum,r,sx,sy,tsx,tsy;

    float **slw, *dip;

    sf_axis az,aa,ad;
    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* input file */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nsx)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dsx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&osx)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&nsy)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dsy)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&osy)) sf_error("No o2= in input");

    if (!sf_histint(in,"n3",&nz)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dz)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&oz)) sf_error("No o3= in input");

    /* output file */
    sf_settype(out,SF_FLOAT);

    if (!sf_getint("na",&na)) na = 360;
    if (!sf_getfloat("da",&da)) da = 0.5;
    if (!sf_getfloat("oa",&oa)) oa = 0.0;

    aa = sf_maxa(na,oa,da);
    sf_oaxa(out,aa,1);

    ad = sf_maxa(1,0,1);
    sf_oaxa(out,ad,2);

    az = sf_maxa(nz,oz,dz);
    sf_oaxa(out,az,3);

    sf_putstring(out,"label1","angle");
    sf_putstring(out,"label2","z");

    if (!sf_getint("nr",&nr)) nr = 2*nsx;
    /* line summation samples */
    if (!sf_getfloat("dr",&dr)) dr = 0.5*dsx;
    /* line summation sampling */
    if (!sf_getfloat("or",&orig)) orig = osx;
    /* line summation origin */

    /* memory allocations */
    slw = sf_floatalloc2(nsx,nsy);
    dip = sf_floatalloc(na);

    for (k = 0; k < nz; k++) {

	sf_floatread(slw[0],nsx*nsy,in);

	for (ia = 0; ia < na; ia++) dip[ia] = 0.0;

	for (ia = 0; ia < na; ia++) {

	    a = (oa + ia*da)*SF_PI/180.0;

	    sum = 0.0;
            nlive = 0;

	    for (i = 0; i < nr; i++) { /* line sum */

		r = orig + i*dr;

		sx = -r*sin(a);
		sy =  r*cos(a);

                /* bilinear interpolation */
                isxm = floorf((sx-osx)/dsx);
                isym = floorf((sy-osy)/dsy);

                if ( (isxm >= 0) && (isym >= 0) && (isxm < (nsx-1)) && (isym < (nsy-1)) ) {

                       tsx = (sx-isxm*dsx-osx)/dsx;
                       tsy = (sy-isym*dsy-osy)/dsy;

                       /* add contribution from this slowness point into line sum */
                       sum += slw[isym][isxm]*(1.-tsx)*(1.-tsy);
                       sum += slw[isym][isxm+1]*tsx*(1.-tsy);
                       sum += slw[isym+1][isxm+1]*tsx*tsy;
                       sum += slw[isym+1][isxm]*(1.-tsx)*tsy;

		       nlive += 1;
                }

	    } /* line sum */

	    if (nlive > 0) dip[ia] = sum/nlive;

	} /* a */

	sf_floatwrite(dip,na,out);

    } /* z */

    exit (0);
}
