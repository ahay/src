/* Objective function of dip estimation with semblance. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
    int i1, i2, n1, n2, np, ip, i, n;
    float p0, dp, p, sum, sum2, f, amp, **xx, *obj;
    sf_file in, of;

    sf_init (argc,argv);
    in = sf_input("in");
    of = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1 in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2 in input");

    if (!sf_getint("np",&np)) np=100;
    /* number of dips */
    if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
    /* dip origin */
    if (!sf_getfloat("dp",&dp)) dp=2*p0/(1.-np);
    /* dip sampling */

    sf_putint(of,"n1",np);
    sf_putfloat(of,"d1",dp);
    sf_putfloat(of,"o1",p0);
    sf_putint(of,"n2",1);

    xx = sf_floatalloc2(n1,n2);
    obj = sf_floatalloc(np);

    sf_floatread (xx[0],n1*n2,in);

    for (ip=0; ip < np; ip++) {
        p = p0 + ip*dp;

	obj[ip]=0.;

	for (i1=0; i1 < n1; i1++) {
	    
	    sum = 0.;
	    sum2 = 0.;
	    n=0;
	    
	    for (i2=0; i2 < n2; i2++) {

		f = i1 + p*i2;
		i = floorf(f);
		
		if (i >=0 && i < n1-1) {
		    f -= i;
		    amp = xx[i2][i]*(1.-f) + xx[i2][i+1]*f;
		    sum += amp;
		    sum2 += amp*amp;
		    n++;
		}
	    }

	    if (n > 0 && sum2 > FLT_EPSILON) obj[ip] += (sum*sum)/(n*sum2);
	}
    }

    sf_floatwrite(obj,np,of);

    exit(0);
}
