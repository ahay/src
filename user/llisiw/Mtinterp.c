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
    int n[3], ns, ns0, nt, is0, s0;
    float o[3], os, os0, d[3], ds, ds0, **t, **tds, *tempt, ss;
    char *type;
    sf_file in, out, deriv, pattern;

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

    if (NULL == (type = sf_getstring("type"))) type="hermit";
    /* type of interpolation (default Hermit) */

    if (type[0] != 'h') {
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
    
    if (NULL != sf_getstring("pattern")) {
	pattern = sf_input("pattern");
    } else {
	pattern = NULL;
    }
    
    if (!sf_getint("ns",&ns0) && 
	(NULL==pattern ||
	 !sf_histint(pattern,"n4",&ns0))) sf_error("Need ns=");
    /* Output source size */
    if (!sf_getfloat("ds",&ds0) && 
	(NULL==pattern ||
	 !sf_histfloat(pattern,"d4",&ds0))) sf_error("Need ds=");
    /* Output source sampling */
    if (!sf_getfloat("os",&os0) &&
	(NULL==pattern ||
	 !sf_histfloat(pattern,"o4",&os0))) sf_error("Need os=");
    /* Output source origin */
    
    sf_putint(out,"n4",ns0);
    sf_putfloat(out,"d4",ds0);
    sf_putfloat(out,"o4",os0);

    /* allocate temporaty memory */
    tempt = sf_floatalloc(nt);

    /* initialization */
    tinterp_init(nt,ds);

    /* loop over sources */
    for (is0=0; is0 < ns0; is0++) {
	s0 = (os0+is0*ds0-os)/ds;
	ss = os0+is0*ds0-os-s0*ds;

	if (s0 < 0) {
	    s0 = 0; ss = 0.;
	}
	if (s0 >= ns-1) {
	    s0 = ns-2; ss = ds;
	}

	/* do interpolation */
	switch (type[0]) {
	    case 'l': /* linear */
		tinterp_linear(tempt,ss,t[s0],t[s0+1]);
		break;

	    case 'p': /* partial */
		tinterp_partial(tempt,ss,n[0],n[1],d[1],t[s0],t[s0+1]);
		break;

	    case 'h': /* hermite */
		tinterp_hermite(tempt,ss,t[s0],t[s0+1],tds[s0],tds[s0+1]);
		break;
	}
	
	sf_floatwrite(tempt,nt,out);
    }

    exit (0);
}
