/* Trapezoid bandpass filter in the frequency domain. 
f1, f2, f3, f4 correspond to four key points of the trapezoid bandpass filter.*/
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

#include <math.h>
#include <rsf.h>
#include "ffilter.h"

int main(int argc, char* argv[])
{
  int  i, n1, n2;
  float *dd,*trace,f[4], d1, fnyq;
  sf_file in, out;
  
  sf_init(argc,argv);
  in  = sf_input("in");
  out = sf_output("out");

  if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input!");
  if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input!");

  fnyq=1/(2*d1);

  if (!sf_getfloat("f1",f))f[0]=0;
  if (!sf_getfloat("f2",f+1))f[1]=5;
  if (!sf_getfloat("f3",f+2))f[2]=fnyq-5;
  if (!sf_getfloat("f4",f+3))f[3]=fnyq;
  
  dd=sf_floatalloc(n1);
  trace=sf_floatalloc(n1);
  n2=sf_leftsize(in,1);
  trapezoid_init(n1,fnyq,f);

  for(i=0;i<n2;i++)
  {
	sf_floatread(dd,n1,in);
 	trapezoid_apply(dd,trace);
        sf_floatwrite(trace,n1,out);
  }

  exit(0);
}



