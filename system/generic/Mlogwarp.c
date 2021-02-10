/* Log warping. */
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
    sf_map4 mo;
    bool inv;
    int n1, i2, n2, i3, n3;
    float o1, d1, o2, d2, eps, t, t0;
    float *trace, *t2, *trace2;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* inversion flag */

    if (inv) {
	if (!sf_histint(in,"n1",&n2)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&d2)) d2=1.;
	if (!sf_histfloat(in,"o1",&o2)) o2=0.;

	if (!sf_histint(in,"n1_logwarp",&n1)) n1=n2;

	if (!sf_getfloat("t0",&t0) && 
	    !sf_histfloat(in,"t0_logwarp",&t0)) sf_error("Need t0=");

	o1 = t0*expf(o2);
	d1 = o2+(n2-1)*d2;
	d1 = (t0*expf(d1)-o1)/(n1-1);
	
	sf_putint(out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);
    } else {
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_getint("pad",&n2)) n2=n1; /* output time samples */

	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	if (!sf_getfloat("t0",&t0)) t0=o1;
	sf_putfloat(out,"t0_logwarp",t0);

	o2 = logf(o1/t0);
	d2 = o1+(n1-1)*d1;
	d2 = (logf(d2/t0) - o2)/(n2-1);

	sf_putint(out,"n1",n2);
	sf_putfloat(out,"d1",d2);
	sf_putfloat(out,"o1",o2);

	sf_putint(out,"n1_t2warp",n1);
    } 

    n3 = sf_leftsize(in,1);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(n2);
    t2 = sf_floatalloc(n2);
    trace2 = sf_floatalloc(n1);

    mo = sf_stretch4_init (n1, o1, d1, n2, eps);

    for (i2=0; i2 < n2; i2++) {
	t = o2+i2*d2;
	t2[i2] = t0*expf(t);
    }    

    sf_stretch4_define (mo,t2,false);
    
    for (i3=0; i3 < n3; i3++) {
	if (inv) {
	    sf_floatread(trace,n2,in);
	    sf_stretch4_apply (false,mo,trace,trace2);
	    sf_floatwrite (trace2,n1,out);
	} else {
	    sf_floatread(trace2,n1,in);
	    sf_stretch4_invert (false,mo,trace,trace2);
	    sf_floatwrite (trace,n2,out);
	}
    }

    exit(0);
}
