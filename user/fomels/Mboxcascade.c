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
    int n1, n2, i1, i2, n, i, m;
    float r, *trace, *trace2;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);

    trace = sf_floatalloc(n1);
    trace2 = sf_floatalloc(n1);
    
    if (!sf_getint("rect",&n)) n=0;
    /* smoothing radius */

    if (!sf_getint("inter",&m)) m=n;
    /* interrupt */

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,inp);
	for (i=1; i < n; i++) {
	    r = -2*cosf(2*SF_PI*i/n);
	    trace2[0] = r*trace[0]+trace[1];
	    for (i1=1; i1 < n1-1; i1++) {
		trace2[i1] = r*trace[i1]+trace[i1-1]+trace[i1+1];
	    }
	    trace2[n1-1] = r*trace[n1-1]+trace[n1-2];
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = (i < n-1)? trace2[i1]: trace2[i1]/(n*n);
	    }
	    if (i==m) break;	
	}
	
	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
