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
#include <rsf.h>

int main (int argc, char* argv[])
{
    int iz,i,j,ia;
    int nsx,nsy,nz,na;
    int nr,nlive,isxm,isym;

    float osx,dsx,osy,dsy,oa,da;
    float or,dr,r,a,line;
    float sx,sy,tsx,tsy;

    float ***slw, **tmp, **dan;

    sf_axis asx,asy,az,aa;
    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    /* input file */
    asx = sf_iaxa(in,1); 
    nsx = sf_n(asx);
    osx = sf_o(asx);
    dsx = sf_d(asx);

    asy = sf_iaxa(in,2); 
    nsy = sf_n(asy);
    osy = sf_o(asy);
    dsy = sf_d(asy);

    az = sf_iaxa(in,3); 
    nz = sf_n(az);

    /* output file */
    sf_settype(out,SF_FLOAT);

    if (!sf_getint("na",&na)) na = 2*floorf(nsx/2) + 1;
    if (!sf_getfloat("da",&da)) da = 90.0/floorf(nsx/2);
    if (!sf_getfloat("oa",&oa)) oa = -90.0;

    aa = sf_maxa(na,oa,da);
    sf_oaxa(out,aa,1);

    sf_oaxa(out,az,2);

    sf_putstring(out,"label1","angle");
    sf_putstring(out,"label2","z");

    /* line summation */
    nr = nsx; 
    or = osx;
    dr = dsx;

    /* memory allocations */
    slw = sf_floatalloc3(nsx,nsy,nz);
    tmp = sf_floatalloc2(nsx,nsy);
    dan = sf_floatalloc2(na,nz);

    sf_floatread(slw[0][0],nsx*nsy*nz,in);

    for (iz = 0; iz < nz; iz++) {

	for (i = 0; i < nsx; i++) { 

	    for (j = 0; j < nsy; j++) {

		tmp[j][i] = slw[iz][j][i];

	    }

	}

	for (ia = 0; ia < na; ia++) {

	    a = (oa + ia*da)*SF_PI/180.0;

	    line = 0.0;
            nlive = 0;

	    dan[iz][ia] = 0.0;

	    for (i = 0; i < nr; i++) {

                /* sum on line */
		r = or + i*dr;

		sx = -r*sin(a);
		sy =  r*cos(a);

                /* bilinear interpolation */
                isxm = floorf((sx-osx)/dsx);
                isym = floorf((sy-osy)/dsy);

                if ( (isxm >= 0) && (isym >= 0) && (isxm < (nsx-1)) && (isym < (nsy-1)) ) {

                       tsx = (sx-isxm*dsx-osx)/dsx;
                       tsy = (sy-isym*dsy-osy)/dsy;

                       /* add contribution from this slowness point into line sum */
                       line += tmp[isym][isxm]*(1.-tsx)*(1.-tsy);
                       line += tmp[isym][isxm+1]*tsx*(1.-tsy);
                       line += tmp[isym+1][isxm+1]*tsx*tsy;
                       line += tmp[isym+1][isxm]*(1.-tsx)*tsy;

		       nlive +=1;
                }

	    } /* sum */

	    if (nlive > 0) dan[iz][ia] = line/nlive;

	} /* a */

    } /* z */

    sf_floatwrite(dan[0],nz*na,out);
	
    exit (0);
}
