/* 2-D irregular data interpolation using natural neighbors and shaping regularization. */
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

#include "nnshape.h"
    
int main(int argc, char* argv[])
{
    bool sym;
    float g1, g2, o3, g3, o1,d1, o2,d2, tol;
    int nd, n1, n2, n12, id, i, three, iter, niter, rect1, rect2, nw;
    float **xy, *z, *m, *m2, *d;
    float xi, dx, dy, dz;
    sf_file in, out, coord, pattern;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    coord = sf_input("coord");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nd)) sf_error("No n1= in input");

    if (NULL != sf_getstring("pattern")) {
	/* pattern file for output dimensions */
	pattern = sf_input("pattern");
	
	if (!sf_histint(pattern,"n1",&n1)) sf_error("No n1= in pattern");
	if (!sf_histint(pattern,"n2",&n2)) sf_error("No n2= in pattern");
	if (!sf_histfloat(pattern,"d1",&d1)) d1=1.;
	if (!sf_histfloat(pattern,"d2",&d2)) d2=1.;
	if (!sf_histfloat(pattern,"o1",&o1)) o1=0.;
	if (!sf_histfloat(pattern,"o2",&o2)) o2=0.;
	
	sf_fileclose(pattern);
    } else {
	if (!sf_getint("n1",&n1)) sf_error("Need n1=");
	if (!sf_getint("n2",&n2)) sf_error("Need n2=");
	if (!sf_getfloat("d1",&d1)) d1=1.;
	if (!sf_getfloat("d2",&d2)) d2=1.;
	if (!sf_getfloat("o1",&o1)) o1=0.;
	if (!sf_getfloat("o2",&o2)) o2=0.;
    }

    n12 = n1*n2;

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"o2",o2);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    xy = sf_floatalloc2(2,nd);
    sf_floatread(xy[0],nd*2,coord);

    if (!sf_getint("rect1",&rect1)) rect1=1;
    if (!sf_getint("rect2",&rect2)) rect2=1;
    /* smoothing regularization */
    
    if (!sf_getint("nw",&nw)) nw=2;
    /* interpolator size */
    
    if (!sf_getbool("sym",&sym)) sym=false;
    /* if y, use symmetric shaping */
    if (!sf_getfloat("tol",&tol)) tol=1e-3;
    /* tolerance for stopping iteration */
    
    nnshape_init(sym,nd, n1,n2, o1,o2, d1,d2, 
		 rect1,rect2, nw, 2, xy);
    sf_gmres_init(n12,niter); 
 
    z = sf_floatalloc (n12);
    m = sf_floatalloc (n12);
    d = sf_floatalloc(nd);

    sf_floatread(d,nd,in);
 
    /* make right-hand side */
    nnshape_back(d,z);
    nnshape_smooth(z);
    /* invert */
    sf_gmres(z,m,nnshape,NULL,niter,tol,true);
    if (sym) nnshape_smooth(m);
    
    sf_floatwrite (m, n12, out);

    exit(0);
}
