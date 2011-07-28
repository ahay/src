/* 2-D regular geometry conversion */
/*
  Copyright (C) 2010 University of Texas at Austin
   
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
#include <math.h>

int main(int argc, char *argv[])
{
    int n1, n2, i2, is, ir, ns, nr, fold;
    float *data, *result, ds, dr, os, r0;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* get data size */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("ns",&ns)) sf_error("Need shot number");
    /* shot number */
    if (!sf_getint("nr",&nr)) sf_error("Need receiver number");
    /* receiver number per shot */

    if (!sf_getfloat("ds",&ds)) ds=0.05;
    /* shot interval */
    if (!sf_getfloat("dr",&dr)) dr=0.025;
    /* receiver interval */

    if (!sf_getfloat("os",&os)) os=0.;
    /* shot origin */
    if (!sf_getfloat("or",&r0)) r0=0.;
    /* receiver origin */

    fold = 1;
    if (ds >= dr) {
	fold = floorf(ds/dr);
    } else {
	sf_error("Need ds >= dr");
    }

    ns = SF_MAX((n1-nr)/2,ns);

    sf_shiftdim(in,out,1);
    sf_putint(out,"n1",nr);
    sf_putint(out,"n2",ns);
    sf_putfloat(out,"d1",dr);
    sf_putfloat(out,"d2",ds);
    sf_putfloat(out,"o1",r0);
    sf_putfloat(out,"o2",os);
    sf_putstring(out,"label1","Offset");
    sf_putstring(out,"label2","Shot");
    
    data = sf_floatalloc(n1);
    result = sf_floatalloc(nr);
    
    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	for (is=0; is < ns; is++) {
	    for (ir=0; ir < nr; ir++) {
		if (is*fold+ir < n1) {
		    result[ir] = data[is*fold+ir];
		} else {
		    result[ir] = 0.;
		}
	    }
	    sf_floatwrite(result,nr,out);
	}
    }

    exit(0);
}
/* 	$Id$	 */
