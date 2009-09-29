/* Local covariance filter */
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
#include <rsf.h>

#include "locov.h"

int main(int argc, char* argv[])
{
    int n1, ia, na, ib, nb;
    float a0, b0, a1, b1, da, db, a, b, *trace, *filt;
    sf_file inp, out;
    
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    
    trace = sf_floatalloc(n1);
    filt = sf_floatalloc(n1);

    if (!sf_getint("na",&na)) na=1;
    if (!sf_getint("nb",&nb)) nb=1;
    if (!sf_getfloat("a0",&a0)) sf_error("Need a0=");
    if (!sf_getfloat("b0",&b0)) sf_error("Need b0=");
    if (!sf_getfloat("a1",&a1)) a1=a0;
    if (!sf_getfloat("b1",&b1)) b1=b0;

    da = (1 >= na)?0.: (a1-a0)/(na-1);
    db = (1 >= nb)?0.: (b1-b0)/(nb-1);

    sf_putint(out,"n2",na);
    sf_putint(out,"n3",nb);

    sf_putfloat(out,"o2",a0);
    sf_putfloat(out,"o3",b0);

    sf_putfloat(out,"d2",da);
    sf_putfloat(out,"d3",db);

    sf_floatread(trace,n1,inp);

    for (ib=0; ib < nb; ib++) {
	b = b0 + ib*db;
	for (ia=0; ia < na; ia++) {
	    a = a0 + ia*da;
	    locov_init(a,b);
	    locov_filt(n1,trace,filt);

	    sf_floatwrite(filt,n1,out);
	}
    }

    exit(0);
}
