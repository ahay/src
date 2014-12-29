/* Add zero trace to original profile in order to improve lateral resolution */

/*
  Copyright (C) 2014 University of Texas at Austin

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
  int  i2, n1, n2, n22, ratio;
  float **dd, d2;
  sf_file in, out;
  
  sf_init(argc,argv);
  in  = sf_input("in");
  out = sf_output("out");

  if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input!");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input!");
  if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input!");
  if (!sf_getint("ratio",&ratio)) ratio=2;

  n22=n2*ratio;
  

  dd=sf_floatalloc2(n1, n22);
  memset(dd[0],0,n1*n22*sizeof(float));

  sf_putint(out,"n2",n22);
  sf_putfloat(out,"d2",d2/ratio);
  for(i2=0;i2<n22;i2+=ratio)
  	sf_floatread(dd[i2],n1,in);

		
  sf_floatwrite(dd[0],n1*n22,out);
  exit(0);
}





