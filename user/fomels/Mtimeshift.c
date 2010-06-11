/* Apply variable time shifts using plane-wave construction. */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

#include "predict.h"

int main(int argc, char* argv[])
{
    int n1, n2, i2, order;
    float *trace, *shift, eps;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    dip = sf_input("dip");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    trace = sf_floatalloc(n1);
    shift = sf_floatalloc(n1);
    predict_init (n1,n2,eps*eps,order,1,false);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,inp);
	sf_floatread(shift,n1,dip);

	predict_step(false,true,trace,shift);

	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
