/* Percentile Data Clipping (C language).
For example: 
A=1,2,3,...,100 
sfpclipc pclip=98 
A'=1,2,3,...,98,98,98	
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
    int i, i3, n1, n2, n3, nthr; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    float t, pclip, *trace=NULL;
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

    if (!sf_getfloat("pclip",&pclip)) pclip=99;
    /* percentile cliping value */

    /* allocate floating point array */
    trace = sf_floatalloc (n1*n2);

    /* loop over gathers */
    for (i3=0;i3<n3;i3++)
    {
        sf_floatread(trace,n1*n2,in);
        nthr=0.5+n1*n2*pclip/100.0;   /* round off A.5=A.0 but A.6!=A.0 */
        t=sf_quantile(nthr-1,n1*n2,trace);

        for (i=0; i < n1*n2; i++)
        {
            if (trace[i] > t) trace[i] = t;
        }
        sf_floatwrite(trace,n1*n2,out);
    }

    exit(0);
}
