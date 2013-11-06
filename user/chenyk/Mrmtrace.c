/* Remove part of traces (resample) in order to make aliasing */

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

int main(int argc, char* argv[])
{
  int  i1, i2, n1, n2, n22, factor;
  float **dd, **ddout, d1, o1, d2, o2, d22, o22;
  sf_file in, out;
  
  sf_init(argc,argv);
  in  = sf_input("in");
  out = sf_output("out");

  if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input!");
  if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input!");
  if (!sf_histfloat(in,"o1",&o1)) o1=0;
  if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input!");
  if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input!");
  if (!sf_histfloat(in,"o2",&o2)) o2=0;

  if (!sf_getint("factor",&factor)) factor=2;
  /* zero part beginning point */
  dd=sf_floatalloc2(n1, n2);
  n22=n2/factor;
  ddout=sf_floatalloc2(n1,n22);

  d22=d2*factor;
  o22=o2+d2*(factor-1);
  sf_floatread(dd[0],n1*n2,in);
  sf_putint(out,"n2",n22);
  sf_putfloat(out,"d2",d22);
  sf_putfloat(out,"o2",o22);

  for(i2=0;i2<n22;i2++)
	for(i1=0;i1<n1;i1++)
	ddout[i2][i1]=dd[(i2+1)*factor-1][i1];

  sf_floatwrite(ddout[0],n1*n22,out);
  exit(0);
}





