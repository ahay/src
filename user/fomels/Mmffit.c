/* Fitting multi-focusing approximations */
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
#include <rsf.h>

#define NC 4

int main(int argc, char* argv[])
{
    const char *type;
    int nh, nm, ih, im, nc;
    float dh, h0, h, dm, m0, m, f, f0, x0, xp, xm, q;
    float fp0, fp, fm0, fm;
    float c[NC], **t, ***dt;
    sf_file in, coef, fit, out;

    sf_init(argc,argv);
    in = sf_input("in");
    coef = sf_input("coef");
    out = sf_output("out");
    fit = sf_output("fit");

    if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&h0)) h0=0.;

    if (!sf_histint(in,"n2",&nm)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dm)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&m0)) m0=0.;

    if (!sf_histint(coef,"n1",&nc) || nc != NC) 
	sf_error("Need n1=%d in coef",NC);
    sf_putint(fit,"n3",NC);
    
    sf_floatread(c,NC,coef);
    sf_fileclose(coef);

    if (!sf_getfloat("x0",&x0)) x0=m0;
    /* central midpoint */

    m0 -= x0;

    if (NULL == (type=sf_getstring("type"))) type="crs";
    /* Type of approximation (crs,mf,nonhyperbolic) */

    t = sf_floatalloc2(nh,nm);
    dt = sf_floatalloc3(nh,nm,NC);

    for (im=0; im < nm; im++) {
	m = m0 + im*dm;

	f0 = c[0] + c[1]*m;
	f = f0*f0 + c[2]*m*m;
	
	for (ih=0; ih < nh; ih++) {
	    h = h0+ih*dh;
	    xp = m+h;
	    xm = m-h;

	    switch (type[0]) {
		case 'c': /* CRS */
		    t[im][ih] = f + c[3]*h*h;
		    dt[0][im][ih] = 2*f0;
		    dt[1][im][ih] = 2*f0*m;
		    dt[2][im][ih] = m*m;
		    dt[3][im][ih] = h*h;
		    break;
		case 'n': /* Nonhyperbolic CRS */
		    fm0 = c[0] + c[1]*xm;
		    fm = fm0*fm0 + c[2]*xm*xm;

		    fp0 = c[0] + c[1]*xp;
		    fp = fp0*fp0 + c[2]*xp*xp;

		    q = sqrtf(fabsf(fp*fm))+SF_EPS;

		    t[im][ih] = 0.5*(f + (2*c[3]+c[1]*c[1]-c[2])*h*h + q);
		    dt[0][im][ih] = f0+0.5*(fp0*fm+fm0*fp)/q;
		    dt[1][im][ih] = f0*m+0.5*(fp0*xp*fm+fm0*xm*fp)/q+c[1]*h*h;
		    dt[2][im][ih] = 0.5*(m*m+0.5*(xp*xp*fm+xm*xm*fp)/q-h*h);
		    dt[3][im][ih] = h*h;
		    break;
		default:
		    sf_error("Unknown type \"%s\"",type); 
		    break;
	    }
	}
    }

    sf_floatwrite(t[0],nh*nm,out);
    sf_floatwrite(dt[0][0],nh*nm*NC,fit);

    exit(0);
}
