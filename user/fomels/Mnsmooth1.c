/* 1-D non-stationary smoothing. */
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

int main (int argc, char* argv[]) 
{
    bool adj;
    int nrep, irep, n1, n2, i1, i2, rect1, shift;
    float* data, **rct, **sft;
    sf_ntriangle tr;
    sf_file in, out, rect, sift;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");
    rect = sf_input("rect");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* Adjoint flag */

    if (NULL != sf_getstring("sift")) {
	sift = sf_input("sift");
    } else {
	sift = NULL;
    }

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    data = sf_floatalloc (n1);

    if (SF_FLOAT != sf_gettype(rect)) sf_error("Need float rect");
    rct = sf_floatalloc2(n1,n2);
    sft = sf_floatalloc2(n1,n2);

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    if (NULL == sift) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		sft[i2][i1] = 0.0f;
	    }
	}
    } else {
	sf_floatread(sft[0],n1*n2,sift);
    }

    sf_floatread(rct[0],n1*n2,rect);

    rect1=1;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    shift = SF_ABS(ceilf(rct[i2][i1]+sft[i2][i1]));
	    if (shift > rect1) rect1=shift;
	}
    }

    tr = sf_ntriangle_init (rect1+1,n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	for (irep=0; irep < nrep; irep++) {
	    if (adj) {
		sf_nsmooth (tr,0,1,false,rct[i2],sft[i2],data);
	    } else {
		sf_nsmooth2 (tr,0,1,false,rct[i2],sft[i2],data);
	    }
	}
	
	sf_floatwrite(data,n1,out);
    }
    
    exit (0);
}

