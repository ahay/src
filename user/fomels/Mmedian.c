/* Compute median on the first axis. */
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

static int n1;
static float *data;

static float middle(float x)
{
    int i1;
    double left, right;

    left = 0.0;
    for (i1=0; i1 < n1/2; i1++ ) {
	left += logf(x-data[i1]);
    }
    right = 0.0;
    for (i1=n1/2; i1 < n1; i1++) {
	right += logf(data[i1]-x);
    }
    return (float) (left-right);
}

int main(int argc, char* argv[]) 
{
    int i2, n2;
    float median, median1, median2, mid, mid1, mid2;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_unshiftdim(in,out,1);

    data = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	if (n1%2) { /* odd number */
	    median = sf_quantile((n1-1)/2,n1,data);
	} else { /* even number */
	    /* Use method from J. Dellinger, 1984, What is the Median of {1,2,3,5} ?: SEP-41, 409-416. */

	    median1 = sf_quantile(n1/2-1,n1,data);
	    median2 = sf_quantile(n1/2,n1,data);

	    if (median1 == median2) {
		median = median1;
	    } else {
                median = 0.5*(median1+median2);
		mid = middle(median);
		mid1 = middle(median1);
		mid2 = middle(median2);

		if (mid > 0.0 && mid1 < 0.0) { 
		    median = sf_zero (middle, median1, median, mid1, mid, SF_EPS, false); 
		} else if (mid < 0.0 && mid2 > 0.0) {
		    median = sf_zero (middle, median, median2, mid, mid2, SF_EPS, false); 
		}
	    }
	}
	sf_floatwrite(&median,1,out);
    }

    exit(0);
}

