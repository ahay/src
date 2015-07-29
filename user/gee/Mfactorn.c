/* Missing data interpolation with 3-D plane-wave filter. */
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

#include "nhelix.h"
#include "nmis2.h"
#include "npfactor.h"

int main(int argc, char* argv[])
{
    bool *mm;
    int ntxy, nt, nx, ny, pt, px, i, niter, npx, npy, miter, n[3], m[3];
    float eps, *pp, *qq, *x;
    nfilter pfilt;
    sf_file in, out, dip, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");
    mask = sf_input("mask");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");
    ntxy = nt*nx*ny;
    n[0] = nt; n[1] = nx; n[2] = ny;

    if (!sf_getfloat("eps",&eps)) eps=0.001;
    if (!sf_getint("nt",&pt)) sf_error("Need nt=");
    if (!sf_getint("nx",&px)) sf_error("Need nx=");
    m[0] = pt; m[1] = px; m[2] = px;

    if (!sf_getint("npx",&npx)) npx=100;
    if (!sf_getint("npy",&npy)) npy=100;

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of factorization iterations */
    if (!sf_getint("miter",&miter)) miter=10;
    /* number of interpolation iterations */

    pp = sf_floatalloc(ntxy);
    qq = sf_floatalloc(ntxy);

    sf_floatread(pp,ntxy,dip);
    sf_floatread(qq,ntxy,dip);
    sf_fileclose(dip);

    pfilt = npfactor(npx, npy, n, m, pp, qq, niter, eps); 
    sf_warning("Done with filters");

    free(pp);
    free(qq);

    x = sf_floatalloc(ntxy);
    mm = sf_boolalloc(ntxy);

    sf_floatread(x,ntxy,mask);

    for (i=0; i < ntxy; i++) {
	mm[i] = (bool) (x[i] != 0.);
    }

    sf_floatread(x,ntxy,in);

    nmis2(niter,ntxy,x,pfilt,mm,true);

    sf_floatwrite(x,ntxy,out);


    exit(0);
}
