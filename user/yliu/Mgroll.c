/* Add linear-chirp ground-roll noise to the data */
/*
  Copyright (C) 2010 University of Texas at Austin
   
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
#include <math.h>

int main(int argc, char *argv[])
{
    int n1, n2, n3, i, j, i3;
    bool rep;
    float *data, d1, d2, o1, o2, max, begf, endf, theta, alpha, beg1, beg2;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* get data size */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");

    n3 = sf_leftsize(in,2);

    if (!sf_getfloat("begf",&begf)) begf=10.;
    /* beginning frequency of ground roll */
    if (!sf_getfloat("endf",&endf)) endf=5.;
    /* ending frequency of ground roll */

    if (!sf_getfloat("theta",&theta)) theta=0.2;
    /* direction of central ground roll */
    if (!sf_getfloat("alpha",&alpha)) alpha=0.1;
    /* width parameter of ground roll */

    if (alpha <= 0) alpha=0.1;
    if (d2 > 0.1) {
	d2 /= 1000.;
	o2 /= 1000.;
    }
    if (d1 > 0.5) {
	d1 /=1000.;
	o1 /=1000.;
    }

    if (!sf_getfloat("beg1",&beg1)) beg1=0.;
    /* radial beginning point at first axis */
    if (!sf_getfloat("beg2",&beg2)) beg2=0.;
    /* radial beginning point at second axis */

    if (!sf_getbool("rep",&rep)) rep=false;
    /* if y, replace data with noise */

    data = sf_floatalloc(n1*n2);

    for (i3=0; i3 < n3; i3++) { 
	sf_floatread(data,n1*n2,in);

	max = 0.;
	for (i=0; i < n1*n2; i++) {
	    if (max <= fabsf(data[i])) max = fabsf(data[i]);
	}
	if (0. == max) max = 1.;
	
	/* loop over traces */
	for (i=0; i < n2; i++) {
	    for (j=0; j < n1; j++) {
		if(rep){
		    data[n1*i+j] =
			max*sinf(2*SF_PI*
				 (begf-(begf-endf)/sqrtf((o1+d1*(n1-beg1))*
							 (o1+d1*(n1-beg1))+
							 (o2+d2*(n2-beg2))*
							 (o2+d2*(n2-beg2)))*
				  sqrtf((o1+d1*(j-beg1))*(o1+d1*(j-beg1))+
					(o2+d2*(i-beg2))*(o2+d2*(i-beg2))))*
				 sqrtf((o1+d1*(j-beg1))*(o1+d1*(j-beg1))+
				       (o2+d2*(i-beg2))*(o2+d2*(i-beg2))))*
			expf(-powf(atanf((o2+d2*(i-beg2))
					/((o1+d1*(j-beg1))+
					  FLT_EPSILON)-theta),2.)/
			     (alpha*alpha));
		} else {
		    data[n1*i+j] += 
			max*sinf(2*SF_PI*
				 (begf-(begf-endf)/sqrtf((o1+d1*(n1-beg1))*
							 (o1+d1*(n1-beg1))+
							 (o2+d2*(n2-beg2))*
							 (o2+d2*(n2-beg2)))*
				  sqrtf((o1+d1*(j-beg1))*(o1+d1*(j-beg1))+
					(o2+d2*(i-beg2))*(o2+d2*(i-beg2))))*
				 sqrtf((o1+d1*(j-beg1))*(o1+d1*(j-beg1))+
				       (o2+d2*(i-beg2))*(o2+d2*(i-beg2))))*
			expf(-powf(atanf((o2+d2*(i-beg2))
					/((o1+d1*(j-beg1))+
					  FLT_EPSILON)-theta),2.)/
			     (alpha*alpha));
		}
	    }
	}
	sf_floatwrite(data,n1*n2,out);
    }

    exit(0);
}
/* 	$Id$	 */
