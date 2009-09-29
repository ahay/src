/* 3-D anisotropic diffusion. */
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
#include "impl3.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, n123, ns;
    float r1, r2, r3, tau, pclip, ***dat=NULL, *dist=NULL;
    bool up, verb;
    char *file=NULL;
    sf_file in=NULL, out=NULL, dst=NULL, snp=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n123 = n1*n2*n3;

    if (!sf_getfloat("rect1",&r1)) sf_error("Need rect1=");
    if (!sf_getfloat("rect2",&r2)) sf_error("Need rect2=");
    if (!sf_getfloat("rect3",&r3)) sf_error("Need rect3=");
    /* smoothing radius */
    if (!sf_getfloat("tau",&tau)) tau=0.1;
    /* smoothing time */
    if (!sf_getfloat("pclip",&pclip)) pclip=50.;
    /* percentage clip for the gradient */
    if (!sf_getbool ("up",&up)) up=false;
    /* smoothing style */
    if (!sf_getbool ("verb",&verb)) verb=false;
    /* verbosity flag */
    if (!sf_getint("nsnap",&ns)) ns=1;
    /* number of snapshots */

    if (NULL != (file = sf_getstring("dist"))) {
	/* inverse distance file (input) */
	dst = sf_input(file);
	dist = sf_floatalloc(n123);
    } else {
	dst = NULL;
	dist = NULL;
    }

    if (NULL != (file = sf_getstring("snap"))) {
	/* snapshot file (output) */
	snp = sf_output(file);
	sf_putint(snp,"n4",ns);
    } else {
	snp = NULL;
    }

    impl3_init (r1, r2, r3, n1, n2, n3, tau, pclip, up, verb, dist, ns, snp);

    dat = sf_floatalloc3(n1,n2,n3);

    if (NULL != dst) sf_floatread(dist,n123,dst);

    sf_floatread (dat[0][0],n123,in);
    impl3_apply (dat,true,false);
    sf_floatwrite (dat[0][0],n123,out);

    exit(0);
}
