/* 2-D implicit finite-difference migration in constant velocity. */
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
#include "fdmig.h"

int main(int argc, char* argv[])
{
    int nx, nz, nw;
    bool hi;
    float dx, dz, dw, vel, beta;
    sf_complex **dat=NULL;
    float **img=NULL;
    sf_file data=NULL, imag=NULL, movie=NULL;

    sf_init(argc, argv);
    data = sf_input("in");
    imag = sf_output("out");

    if (SF_COMPLEX != sf_gettype(data)) sf_error("Need complex input");
    if (!sf_histint(data,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(data,"n2",&nw)) sf_error("No n2= in input");
    if (!sf_histfloat(data,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(data,"d2",&dw)) sf_error("No d2= in input");

    if (!sf_getint("nz",&nz)) nz=2*(nw-1);
    /* vertical time samples */
    if (!sf_getfloat("dz",&dz)) dz=1./(nz*dw);
    /* vertical time sampling */

    sf_settype(imag,SF_FLOAT);
    sf_putint(imag,"n2",nz);
    sf_putfloat(imag,"d2",dz);

    if (NULL != sf_getstring("movie")) {
	movie = sf_output("movie");
	sf_putint(movie,"n2",nz);
	sf_putfloat(movie,"d2",dz);
	sf_putint(movie,"n3",nw);
	sf_putfloat(movie,"d3",dw);
    } else {
	movie = NULL;
    }

    if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");
    /* constant velocity */

    if (!sf_getbool("hi",&hi)) hi=true;
    /* if y, use 45-degree; n, 15-degree */

    if (!sf_getfloat("sixth",&beta)) beta=1./12;
    /* one-sixth trick */

    fdmig_init(hi, nx, nz, nw, dx, dz, dw, vel, beta);

    dat = sf_complexalloc2(nx,nw);
    img = sf_floatalloc2(nx,nz);

    sf_complexread(dat[0],nx*nw,data);
    fdmig (dat, img, movie);
    sf_floatwrite(img[0],nx*nz,imag);

    exit(0);
}
