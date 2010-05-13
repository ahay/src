/* Recursive Gaussian smoothing on the fast axis. */
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

#include "recgauss.h"

int main (int argc, char* argv[]) 
{
    bool der;
    int i1, n1, i2, n2, irep, nrep;
    float *data, *data2, rect;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getbool("der",&der)) der=false;
    /* compute derivative */

    if (!sf_getfloat("rect",&rect)) rect=1;
    /* smoothing radius */ 
    recgauss_init (n1,der,rect);

    data = sf_floatalloc (n1);
    data2 = der? sf_floatalloc (n1): NULL;

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    for (i2=0; i2 < n2; i2++) {
	if (der) {
	    sf_floatread(data2,n1,in);
	    data[0] = data2[1]-data2[0];
	    for (i1=1; i1 < n1-1; i1++) {
		data[i1] = 0.5*(data2[i1+1]-data2[i1-1]);
	    }
	    data[n1-1] = data2[n1-1]-data2[n1-2];
	} else {
	    sf_floatread(data,n1,in);
	}

	for (irep=0; irep < nrep; irep++) {
	    recgauss (data);
	}
	sf_floatwrite(data,n1,out);
    }    

    exit (0);
}


