/* Zero part of traces in order to make aliasing */

/*
  Copyright (C) 2013 Yangkang Chen, University of Texas at Austin

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
  int  i1, i2, i3, ijump, ilast, njump, jump, last, n1, n2, n3;
  float **dd, d1, o1, d2, o2, beg;
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

  if (!sf_getfloat("beg",&beg)) beg=o2;
  /* zero part beginning point */
  if (!sf_getint("j", &jump)) jump=2;
  /* jump step between two consecutive zero parts */
  if (!sf_getint("l", &last)) last=1;
  /* length of each zero part */

  njump=(o2+(n2-1)*d2)/jump+1;

  dd=sf_floatalloc2(n1, n2);

  n3 = sf_leftsize(in,2); 
  for(i3=0;i3<n3;i3++)
{  sf_floatread(dd[0],n1*n2,in);


  for(ijump=0;ijump<njump;ijump++)
     for(ilast=0;ilast<last;ilast++)
	 { i2=beg+ijump*jump+ilast;
	     for(i1=0;i1<n1;i1++)
		dd[i2][i1]=0;}
  sf_floatwrite(dd[0],n1*n2,out);}
  exit(0);
}





