/* Kirchhoff 2-D post-stack least-squares time migration with sparse constrains. 

 Antialiasing by reparameterization. */
/*
  Copyright (C) 2013 China University of Petroleum (East China)
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

#include <math.h>

#include <rsf.h>

#include "kirchnew.h"
#include "inverts.h"

int main(int argc, char* argv[])
{
    int n12, n1, n2, n3, i3, sw, niter, liter;
    bool hd, ps, verb;
    float *data, *modl, *vrms, *error=NULL;
    float o1,d1,o2,d2,eps;
    char *errfile;
    sf_file in, out, vel, err=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("hd",&hd)) hd=true;
    /* if y, apply half-derivative filter */
    if (!sf_getbool("ps",&ps)) ps=true;
    /* if y, apply pseudo-unitary weighting */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    if (!sf_getint("sw",&sw)) sw=0;
    /* if > 0, select a branch of the antialiasing operation */
    if (!sf_getint("niter",&niter)) niter=5;
    /* number of non-linear iterations, when niter=1, it's linear */
    if (!sf_getint("liter",&liter)) liter=5;
    /* number of linear iterations */
    if (!sf_getfloat("eps",&eps)) eps=0.;
    /* regularization parameters */

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

    vrms = sf_floatalloc(n1);

    vel = sf_input("velocity");
    sf_floatread(vrms,n1,vel);
    sf_fileclose(vel);
    
    n12 = n1*n2;
    data = sf_floatalloc(n12);
    modl = sf_floatalloc(n12);

    if (niter > 0) {
	errfile = sf_getstring("err");
	/* output file for error */
	if (NULL != errfile) {
	    err = sf_output(errfile);
	    sf_putint(err,"n1",liter);
	    sf_putfloat(err,"d1",1);
	    sf_putfloat(err,"o1",1);
	    sf_putstring(err,"label1","Iteration Number");
	    sf_putstring(err,"label2","Relative Squared Error");
	    sf_putint(err,"n2",1);
	}
	error = sf_floatalloc(liter);
    }

    kirchnew_init (vrms, o1, d1, d2, n1, n2, sw, ps, hd);

    for (i3=0; i3 < n3; i3++) {

	sf_floatread (data,n12,in);

	inverts(kirchnew_lop,niter,niter,liter,n12,n12,modl,data,error,verb,eps);

	sf_floatwrite (modl,n12,out);

	if (NULL != err) sf_floatwrite(error,liter,err);
    }

    sf_warning(".");

    exit(0);
}
