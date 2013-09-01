/* Image-domain waveform tomography (image perturbation). */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
    int i1, n1, i2, n2, ih, nh, idx;
    float *slope, *deriv, thres;
    float *idz, *idh;
    sf_file in, out, pz, ph;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getfloat("thres",&thres)) thres=0.01;
    /* slope thresholding */

    /* read input */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");
    if (!sf_histint(in,"n3",&nh)) sf_error("No n3= in input.");

    slope = sf_floatalloc(n1*n2*nh);
    sf_floatread(slope,n1*n2*nh,in);
    
    /* read image derivatives */
    if (NULL == sf_getstring("pz"))
	sf_error("Need dI/dz pz=");
    pz = sf_input("pz");

    idz = sf_floatalloc(n1*n2*nh);
    sf_floatread(idz,n1*n2*nh,pz);
    sf_fileclose(pz);

    if (NULL == sf_getstring("ph"))
	sf_error("Need dI/dh ph=");
    ph = sf_input("ph");

    idh = sf_floatalloc(n1*n2*nh);
    sf_floatread(idh,n1*n2*nh,ph);
    sf_fileclose(ph);

    /* allocate temporary memory */
    deriv = sf_floatalloc(n1*n2*nh);

    /* loop over image point */
    for (ih=0; ih < nh; ih++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		idx = ih*n1*n2+i2*n1+i1;
		
		if (fabsf(slope[idx]) < thres || ih == (nh-1)/2) {
		    deriv[idx] = 0.;
		    continue;
		}		

		deriv[idx] = (idz[idx]-slope[idx]*idh[idx])
		    /sqrtf(1.+slope[idx]*slope[idx]);

		if (slope[idx] > 0. && ih > (nh-1)/2) 
		    deriv[idx] = -deriv[idx];
		if (slope[idx] < 0. && ih < (nh-1)/2) 
		    deriv[idx] = -deriv[idx];
	    }
	}
    }

    /* write output */
    sf_floatwrite(deriv,n1*n2*nh,out);

    exit(0);
}
