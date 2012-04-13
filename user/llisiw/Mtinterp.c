/* Traveltime interpolation by cubic Hermite spline */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include "tinterp.h"

int main (int argc,char* argv[]) 
{
    int n[3], ns, ns0, nt, interp, is0, ss, s0;
    float o[3], os, d[3], ds, ds0, **t, **tds, *tempt;
    char *what;
    sf_file in, out, deriv;

    sf_init (argc, argv);
    in  = sf_input("in");
    out = sf_output("out");

    /* read input dimensions */
    if(!sf_histint(in,"n1",n  )) sf_error("No n1= in input");
    if(!sf_histint(in,"n2",n+1)) sf_error("No n2= in input");
    if(!sf_histint(in,"n3",n+2)) n[2]=1;
    if(!sf_histint(in,"n4",&ns)) sf_error("No ns= in input");

    if (ns <= 1) sf_error("Provide at least two shots");

    if(!sf_histfloat(in,"d1",d  )) sf_error("No d1= in input");
    if(!sf_histfloat(in,"d2",d+1)) sf_error("No d2= in input");
    if(!sf_histfloat(in,"d3",d+2)) d[2]=d[1];
    if(!sf_histfloat(in,"d4",&ds)) ds=d[1];

    if(!sf_histfloat(in,"o1",o  )) o[0]=0.;
    if(!sf_histfloat(in,"o2",o+1)) o[1]=0.;
    if(!sf_histfloat(in,"o3",o+2)) o[2]=o[1];
    if(!sf_histfloat(in,"o4",&os)) os=o[1];

    /* read traveltime file */
    nt = n[0]*n[1]*n[2];
    
    t = sf_floatalloc2(nt,ns);
    sf_floatread(t[0],nt*ns,in);

    if (NULL == (what = sf_getstring("what"))) what="expanded";
    /* Hermite basis functions (default expanded) */

    if (what[0] == 'l') {
	/* linear interpolation */
	deriv = NULL;
	tds = NULL;
    } else {
	/* read derivative file */
	if (NULL == sf_getstring("deriv"))
	    sf_error("Need derivative deriv=");
	deriv = sf_input("deriv");

	tds = sf_floatalloc2(nt,ns);
	sf_floatread(tds[0],nt*ns,deriv);
	sf_fileclose(deriv);
    }
    
    if (!sf_getint("interp",&interp)) interp=1;
    /* number of interpolation */
    
    ds0 = ds/(float)(1+interp);
    ns0 = (ns-1)*interp+ns;

    sf_putint(out,"n4",ns0);
    sf_putfloat(out,"d4",ds0);

    /* allocate temporaty memory */
    tempt = sf_floatalloc(nt);

    /* initialization */
    tinterp_init(nt,ds,what);

    /* loop over sources */
    for (is0=0; is0 < ns0; is0++) {
	ss = is0%(interp+1);
	s0 = (is0-ss)/(interp+1);

	/* continue if input sources */
	if (ss == 0) {
	    sf_floatwrite(t[s0],nt,out);
	    continue;
	}
	
	/* do interpolation */
	switch (what[0]) {
	    case 'l': /* linear */
		tinterp_linear(tempt,((float)ss)*ds0,t[s0],t[s0+1]);
		break;

	    default: /* cubic Hermite spline */
		tinterp_hermite(tempt,((float)ss)*ds0,t[s0],t[s0+1],tds[s0],tds[s0+1]);
		break;
	}
	
	sf_floatwrite(tempt,nt,out);
    }

    exit (0);
}
