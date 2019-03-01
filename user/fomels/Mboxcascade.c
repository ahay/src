/* Box filter cascade */
/*
  Copyright (C) 2019 University of Texas at Austin
  
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
    bool tri;
    int n1, n2, i1, i2, n, i, m;
    float r;
    sf_complex *trace, c, new, pre;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);

    trace = sf_complexalloc(n1);

    if (!sf_getint("rect",&n)) n=1;
    /* smoothing radius */

    if (!sf_getint("inter",&m)) m=n;
    /* interrupt */

    if (!sf_getbool("triangle",&tri)) tri=false;
    /* if triangle filter */

    for (i2=0; i2 < n2; i2++) {
	sf_complexread(trace,n1,inp);
	for (i=1; i < n; i++) {
	    r = 2*SF_PI*i/n;
	    c = -1.0f/sf_cmplx(cosf(r),sinf(r));
	    pre = 0.0f;
	    for (i1=0; i1 < n1; i1++) {
		new = trace[i1];
		trace[i1] += c*pre;
		pre = new;
		if (i==n-1) trace[i1] /= n;
	    }
	    if (tri) {
		pre = 0.0f;
		for (i1=n1-1; i1 >= 0; i1--) {
		    new = trace[i1];
		    trace[i1] += c*pre;
		    pre = new;
		    if (i==n-1) trace[i1] /= n;
		}
	    }
	    if (i==m) break;	
	}
	
	sf_complexwrite(trace,n1,out);
    }

    exit(0);
}
