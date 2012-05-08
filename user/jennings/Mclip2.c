/* One- or two-sided data clipping.

sfclip2 is a generalization of sfclip.

Clip values above xu:         sfclip2 < in.rsf > out.rsf upper=xu
Clip values below xl:         sfclip2 < in.rsf > out.rsf lower=xl
Clip values outside [xu,xl]:  sfclip2 < in.rsf > out.rsf upper=xu lower=xl

sfclip2 < in.rsf > out.rsf upper=x lower=-x

is equivalent to

sfclip < in.rsf > out.rsf clip=x
*/

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
#include <float.h>

int main(int argc, char* argv[])
{
    off_t i, n, nbuf;
    float upper, lower, *trace=NULL;
    sf_file in=NULL, out=NULL; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    n = sf_filesize(in);

    if (!sf_getfloat("upper",&upper)) upper=+FLT_MAX;
    /* upper clip value */

    if (!sf_getfloat("lower",&lower)) lower=-FLT_MAX;
    /* lower clip value */

    /* allocate floating point array */
    nbuf = BUFSIZ/sizeof(float);
    trace = sf_floatalloc (nbuf);

    /* loop over traces */
    for (n = sf_filesize(in); n > 0; n -= nbuf)
    {
        if (nbuf > n) nbuf=n;

        sf_floatread(trace,nbuf,in);

        for (i=0; i < nbuf; i++)
        {
            if (trace[i] > upper) trace[i] = upper;
            if (trace[i] < lower) trace[i] = lower;
        }

        sf_floatwrite(trace,nbuf,out);
    }



    exit(0);
}
