/* Circular statistics correlation for 2D data. */
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

int main(int argc, char* argv[])
{

    int nbin, iax, nx, ny;
    float mean, strength, sigma;
    float dx, dy, xi, yi;
    float **data; 
    float *circ;
    bool verb;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ny)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&xi)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d2",&dy)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&yi)) sf_error("No o2= in input");

    if (!sf_getint("nbin",&nbin)) nbin=10;
    /* number of distance bins */
    if (!sf_getint("iax",&iax)) iax=1;
    /* flag for data type 1=angular 2=axial */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    data = sf_floatalloc2(nx,ny);
    circ = sf_floatalloc(nbin);

    sf_floatread(data[0],nx*ny,in);

    if (iax = 2) 
    {
     for (i = 0; i < ndata; i++) axes[i] *= 2.0;
    }

    momentcircular(nx,ny,data,mean,strength);
    sigma = 1.0 - strength*strength;
    if (iax = 2) mean *= 0.5;

    if (verb) sf_warning("got mean=%g strength=%g sigma=%g",mean*180.*pi,strength,sigma);

    binning(data,nbin,nx,ny,dx,dy,xi,yi,strength,verb);

    sf_floatwrite (circ,nbin,out);

    exit (0);
}

/* 	$Id$	 */

