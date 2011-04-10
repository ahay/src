/* Offset continuation by log-stretch F-K operator. */
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
    int nw, nk, iw, ik, nh, ih;
    sf_complex *oper=NULL;
    float dw,dk, ow,ok, k,w,h, h0, dh, eps,amp,phase;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nk)) sf_error("No n2= in input");

    if (!sf_histfloat (in,"d1",&dw)) sf_error("No d1= in input");
    if (!sf_histfloat (in,"o1",&ow)) sf_error("No o1= in input");
    if (!sf_histfloat (in,"d2",&dk)) sf_error("No d2= in input");
    if (!sf_histfloat (in,"o2",&ok)) sf_error("No o2= in input");

    if (!sf_getfloat("h",&h)) sf_error("Need h=");
    /* final offset */
    if (!sf_getint("nh",&nh)) nh=1;
    /* number of offset steps */
    if (!sf_getfloat("h0",&h0)) h0=0.;
    /* initial offset */
    dh = (h-h0)/nh;

    sf_putint(out,"n3",nh);
    sf_putfloat(out,"d3",dh*2.);
    sf_putfloat(out,"o3",(h0+dh)*2.);

    oper = sf_complexalloc (nw);

    for (ih = 0; ih < nh; ih++) {
	h = h0 + (ih+1)*dh;

	for (ik = 0; ik < nk; ik++) {
	    k = ok+ik*dk;

	    for (iw = 0; iw < nw; iw++) {
		w = ow + iw*dw;

		if (fabsf (w) > FLT_EPSILON) {
		    eps = 2.*k*h/w;
		    eps = sqrtf (1+eps*eps);
		    amp = sqrtf(0.5*(1/eps+1.))*expf(0.5*(1-eps));
		    phase = 1-eps+logf(0.5*(1+eps));
		    phase *= -SF_PI*w;

		    oper[iw] = sf_cmplx(amp*cosf(phase),amp*sinf(phase));
		} else {
		    oper[iw] = sf_cmplx(0.,0.);
		}
	    }

	    sf_complexwrite(oper,nw,out);
	}
    }

    exit(0);
}
