/* Soft clip the data.
Uses softplus function
Performs lower clipping then upper clipping if both specified

lower clipping:
y = output, x = input, k = sharpness, c = lower clip value
y = ln(1+exp(k*(x-c)))/k + c
 
upper clipping:
y = output, x = input, k = sharpness, c = upper clip value
y = -ln(1+exp(k*(c-x)))/k + c
*/
/*
  Copyright (C) 2020 University of Texas at Austin

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
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int i, n, nbuf;
    float upper, lower, sharp, *trace;

    sf_file in=NULL, out=NULL; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
	
	bool clipping_upper = true ;

    if (!sf_getfloat("upper",&upper)) {upper=0.0; clipping_upper = false;}
    /* upper clip value */
	
	bool clipping_lower = true;
    if (!sf_getfloat("lower",&lower)) {lower=0.0; clipping_lower=false;}
    /* lower clip value */
	
    if (!sf_getfloat("sharp",&sharp)) sharp = 1.0/fmax(upper*upper,lower*lower);
	
	if (sharp != sharp){ sharp = 1.0;}
    /* sharpness */	
	if( sharp <= 0.0 ){
		sf_error("Sharpness must be > 0 ");
	}

    /* allocate floating point buffer */
    nbuf = sf_bufsiz(in)/sizeof(float);
    trace = sf_floatalloc (nbuf);

    /* process buffers */
    for (n = sf_filesize(in); n > 0; n -= nbuf) {
	if (nbuf > n) nbuf=n;

	sf_floatread(trace,nbuf,in);

	for (i=0; i < nbuf; i++) {
		if (clipping_lower){
			trace[i] = log(1.0+exp((trace[i]-lower)*sharp)) / sharp + lower;
		}
		if (clipping_upper){
			trace[i] = -1*log(1.0+exp((upper-trace[i])*sharp))/sharp + upper;
		}
	}

	sf_floatwrite(trace,nbuf,out);
    }


    exit(0);
}
