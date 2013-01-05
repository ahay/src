/* 2D divn */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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
#include "runtime.h"

int main(int argc, char* argv[])
{
    bool verb;
    int i1, i3, n12, n3, niter, n[2], rect[2];
    float *u1, *u2, *u3, norm;
    sf_file inp, out, den;

    sf_init(argc,argv);
    inp = sf_input("in");
    den = sf_input("den");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    
    if (!sf_histint(inp,"n1",n)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",n+1)) sf_error("No n2= in input");
    n3 = sf_leftsize(inp,2);
    n12 = n[0]*n[1];

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* smoothing radius */

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */

    u1 = sf_floatalloc(n12);
    u2 = sf_floatalloc(n12);
    u3 = sf_floatalloc(n12);

	runtime_init(n12*sizeof(float));
    for (i3=0; i3 < n3; i3++) 
	{
		sf_floatread(u1, n12, inp);
		sf_floatread(u2, n12, den);

    	sf_divn_init(2, n12, n, rect, niter, verb);
		/* smooth division */
		norm = 0.;
		for (i1=0; i1 < n12; i1++) 
		    norm += u2[i1] * u2[i1];
		norm = sqrtf(n12/norm);

		for (i1=0; i1 < n12; i1++) 
		{
			u1[i1] *= norm;
			u2[i1] *= norm;
		}
		sf_divn (u1, u2, u3);

		sf_floatwrite(u3, n12, out);

		sf_divn_close();
		norm = runtime(1);
		sf_warning("%d of %d, %f MB/sec;", i3, n3, norm);
	}

	free(u1);
	free(u2);
	free(u3);
    exit(0);
}


