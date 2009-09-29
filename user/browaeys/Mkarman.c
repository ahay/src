/* Estimation of von Karman autocorrelation 1D spectrum. */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "karman.h"

int main(int argc, char* argv[])
{
    int nk, i2, n2, niter, nliter;
    float k0, dk, x0, x[3], *data, *pred;

    bool verb;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nk)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&dk)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&k0)) sf_error("No o1= in input");

    if (!sf_getfloat("x0",&x0)) x0=1.;
    /* initial squared length scale */
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getint("nliter",&nliter)) nliter=1;
    /* number of reweighting iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    data = sf_floatalloc(nk);
    pred = sf_floatalloc(nk);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,nk,in);

	x[0] = x0;
	x[1] = x[2] = 1.;
	
	karman(data,pred,nliter,niter,x,nk,dk,k0,verb);
        
	sf_floatwrite (pred,nk,out);
    }

    exit (0);
}

/* 	$Id$	 */

