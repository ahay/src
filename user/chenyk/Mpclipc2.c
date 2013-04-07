/* One-or Two-sided Percentile Data clipping (C language).
For example: 
A=1,2,3,...,100 
sfpclipc2 upclip=98 lpclip=3 
A'=3,3,3,...,98,98,98	
*/

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
#include <float.h>

int main(int argc, char* argv[])
{
    int i, i3, n1, n2, n3, nthr1, nthr2; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    float tu, tl, upclip, lpclip, *trace=NULL;
    sf_file in=NULL, out=NULL; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if(!sf_histint(in,"n1",&n1)) sf_error("No n1 in input!");
    if(!sf_histint(in,"n2",&n2)) sf_error("No n2 in input!");

    n3=sf_leftsize(in,2);

    if (!sf_getfloat("upclip",&upclip)) upclip=99;
    /* percentile upper cliping value */
    if (!sf_getfloat("lpclip",&lpclip)) lpclip=0;
    /* percentile lower cliping value */

    /* allocate floating point array */
    trace = sf_floatalloc (n1*n2);

    /* loop over gathers */
    for (i3=0;i3<n3;i3++)
    {
        sf_floatread(trace,n1*n2,in);
        nthr1=0.5+n1*n2*upclip/100.0;
        nthr2=0.5+n1*n2*lpclip/100.0;

        tu=sf_quantile(nthr1-1,n1*n2,trace);
        tl=sf_quantile(nthr2-1,n1*n2,trace);

        for (i=0; i < n1*n2; i++)
        {
            if (trace[i] > tu) trace[i] = tu;
	    if (trace[i] < tl) trace[i] = tl;
        }
        sf_floatwrite(trace,n1*n2,out);
    }

    exit(0);
}
