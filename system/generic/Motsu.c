/* Compute a threshold value from histogram using Otsu's algorithm. */
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
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[]) 
{
    int i1, n1, *hist, thr, w1, w2, s1, sum, total;
    float var, varmax, m1, m2, o1, d1;
    sf_file in;

    sf_init(argc,argv);
    in = sf_input("in");

    if (SF_INT != sf_gettype(in)) sf_error("Need int input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.0f;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.0f;
    
    hist = sf_intalloc(n1);
    sf_intread(hist,n1,in);

    total = sum = 0;
    for (i1=0; i1 < n1; i1++) {
	total += hist[i1];
	sum   += i1*hist[i1];
    }

    w1 = w2 = s1 = thr = 0;
    varmax = 0.0f;
    for (i1=0; i1 < n1; i1++) {
	w1 += hist[i1];
	if (0 == w1) continue;

	w2 = total-w1;
	if (0 == w2) break;
	
	s1 += i1*hist[i1];

	/* mean values */
	m1 = (float) s1/w1;         
	m2 = (float) (sum - s1)/w2; 

	/* between class variance */
	var = w1*w2*(m1-m2)*(m1-m2);

	if (var > varmax) {
	    varmax = var;
	    thr = i1;
	}
    }

    printf("threshold=%g\n\n",o1+(thr+0.5)*d1);
    exit(0);
}
