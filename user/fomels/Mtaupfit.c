/* Fitting tau-p approximations */
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

int main(int argc, char* argv[])
{
    const char *type;
    int np, ip, nc, i2, n2, nc2;
    float dp, p0, p;
    float *c, *t, **dt;
    sf_file in, coef, fit, out;

    sf_init(argc,argv);
    in = sf_input("in");
    coef = sf_input("coef");
    out = sf_output("out");
    fit = sf_output("fit");

    if (!sf_histint(in,"n1",&np)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dp)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&p0)) p0=0.;
    n2 = sf_leftsize(in,1);

    if (NULL == (type=sf_getstring("type"))) type="iso";
    /* Type of approximation (iso,vti) */

    nc2 = (type[0]=='i')? 2:3;

    if (!sf_histint(coef,"n1",&nc) || nc != nc2) 
	sf_error("Need n1=%d in coef",nc2);

    c = sf_floatalloc(nc);

    sf_putint(fit,"n2",nc);
    sf_shiftdim(in, fit, 2);

    t = sf_floatalloc(np);
    dt = sf_floatalloc2(np,nc);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(c,nc,coef);

	for (ip=0; ip < np; ip++) {
	    p = p0 + ip*dp;
	    p *= p;
	    
	    switch (type[0]) {
		case 'i': 
		    dt[0][ip] = 1;
		    dt[1][ip] = -p;
		    t[ip] = c[0]-p*c[1];
		    break;
		case 'v':
		    dt[0][ip] = (1-p*c[1])/(1-p*c[2]);
		    t[ip] = c[0]*dt[0][ip];
		    dt[1][ip] = -p*c[0]/(1-p*c[2]);
		    dt[2][ip] = p*t[ip]/(1-p*c[2]);
		    break;
		default:
		    sf_error("Unknown type \"%s\"",type); 
		    break;
	    }
	}

	sf_floatwrite(t,np,out);
	sf_floatwrite(dt[0],np*nc,fit);
    }

    exit(0);
}
