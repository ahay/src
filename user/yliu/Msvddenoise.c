/* SVD denoising */
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

#include "svddenoise.h"

int main(int argc, char *argv[])
{
    int i, n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    float *pp, *qq;
    sf_file in, out;
    float pclip;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/

    if (!sf_getfloat("pclip",&pclip)) pclip=99.;
    /* data clip percentile (default is 99)*/
    if (pclip <=0. || pclip > 100.)
	sf_error("pclip=%g should be > 0 and <= 100",pclip);
 
    pp = sf_floatalloc(n1*n2);
    qq = sf_floatalloc(n1*n2);

    for(i=0;i<n3;i++)  {
	sf_floatread(pp,n1*n2,in);
	svd_denoise(n2,n1,pclip,pp,qq);
	sf_floatwrite(qq,n1*n2,out);
    }

    exit(0);
}

/* 	$Id$	 */
