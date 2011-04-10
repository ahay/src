/* Computes Azimuth Move-Out (AMO) operator in the f-k log-stretch domain */
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

#include <math.h>
#include <float.h>
#include <rsf.h>

int main(int argc, char* argv[])
{
    int nw, nx, ny, iw, ix, iy;
    sf_complex *oper=NULL;
    float dw,dx,dy, ow,ox,oy, w,x,y, x1,x2, h1,h2,f1,f2, maxe;
    float eps1,eps2,amp1,amp2,phase1,phase2,amp;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");

    if (!sf_histfloat (in,"o1",&ow)) sf_error("No o1= in input");
    if (!sf_histfloat (in,"d1",&dw)) sf_error("No d1= in input");
    if (!sf_histfloat (in,"o2",&ox)) sf_error("No o2= in input");
    if (!sf_histfloat (in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat (in,"o3",&oy)) sf_error("No o3= in input");
    if (!sf_histfloat (in,"d3",&dy)) sf_error("No d3= in input");

    if (!sf_getfloat("h1",&h1)) sf_error("Need h1=");
    /* input offset */
    if (!sf_getfloat("h2",&h2)) sf_error("Need h2=");
    /* output offset */
    if (!sf_getfloat("f1",&f1)) sf_error("Need f1=");
    /* input azimuth in degrees */
    if (!sf_getfloat("f2",&f2)) sf_error("Need f2=");
    /* output azimuth in degrees */

    if (!sf_getfloat("maxe",&maxe)) maxe=10.;
    /* stability constraint */

    f1 *= SF_PI/180.;
    f2 *= SF_PI/180.;

    oper = sf_complexalloc (nw);

    for (iy=0; iy < ny; iy++) {
	y = oy + iy*dy;
	for (ix=0; ix < nx; ix++) {
	    x = ox + ix*dx;
	    x1 = x*cosf(f1) + y*sinf(f1);
	    x2 = x*cosf(f2) + y*sinf(f2);
	    for (iw=0; iw < nw; iw++) {
		w = ow + iw*dw;
		if (fabsf (w) > FLT_EPSILON) {
		    eps1 = 2.*fabsf(x1*h1/w);
		    eps2 = 2.*fabsf(x2*h2/w);
		    if (eps1 <= maxe && eps2 <= maxe) {
			eps1 = hypotf (1.,eps1);
			eps2 = hypotf (1.,eps2);

			amp1 = 1./eps1+eps1;
			amp2 = 1./eps2+eps2;
			phase1 = 1-eps1+logf(0.5*(1.+eps1));
			phase2 = 1-eps2+logf(0.5*(1.+eps2));

			amp = expf(0.5*(eps1-logf(amp1)+logf(amp2)-eps2));

			phase1 = SF_PI*(phase2-phase1)*w;
			oper[iw] = sf_cmplx(amp*cosf(phase1),amp*sinf(phase1));
		    } else {
			oper[iw] = sf_cmplx(0.,0.);
		    }
		} else {
		    oper[iw] = sf_cmplx(0.,0.);
		}
	    }
	    sf_complexwrite (oper,nw,out);
	}
    }

    exit (0);
}
