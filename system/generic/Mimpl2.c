/* 2-D anisotropic diffusion. */
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
    int n1, n2, n12, n3, i3, ns;
    float r1, r2, tau, pclip, **dat=NULL, *dist=NULL, **dat2=NULL;
    bool up, verb, lin, adj;
    char *file=NULL;
    sf_file in=NULL, out=NULL, dst=NULL, snp=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_getfloat("rect1",&r1)) sf_error("Need rect1=");
    /* vertical smoothing */
    if (!sf_getfloat("rect2",&r2)) sf_error("Need rect2=");
    /* horizontal smoothing */
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
    if (!sf_getbool ("lin",&lin)) lin=false;
    /* if linear operator */
    if (!sf_getbool ("adj",&adj)) adj=false;
    /* adjoint flag */

    if (NULL != (file = sf_getstring("dist"))) {
	/* inverse distance file (input) */
	dst = sf_input(file);
	dist = sf_floatalloc(n12);
    } else {
	dst = NULL;
	dist = NULL;
    }

    if (NULL != (file = sf_getstring("snap"))) {
	/* snapshot file (output) */
	snp = sf_output(file);
	sf_putint(snp,"n3",ns);
    } else {
	snp = NULL;
    }

    sf_impl2_init (r1, r2, n1, n2, tau, pclip, up, verb, dist, ns, snp);

    dat = sf_floatalloc2(n1,n2);
    dat2 = lin? sf_floatalloc2(n1,n2): NULL;

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d;",i3,n3);
	if (NULL != dst) sf_floatread(dist,n12,dst);

	sf_floatread (dat[0],n12,in);

	if (lin) {
	    if (adj) {
		sf_impl2_lop(true,false,n12,n12,dat2[0],dat[0]);
	    } else {
		sf_impl2_lop(false,false,n12,n12,dat[0],dat2[0]);
	    }
	    sf_floatwrite (dat2[0],n12,out);
	} else {
	    sf_impl2_apply (dat,true,false);
	    sf_floatwrite (dat[0],n12,out);
	}
    }
    sf_warning(".");

    exit(0);
}
