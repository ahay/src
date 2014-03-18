/* Matlab-like matrix fliplr/flipud (flip left and right or up and down)
*/
/*
  Copyright (C) 2014 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
   
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

int main(int argc, char *argv[])
{
    char *mode;
    int n1, n2, i1, i2;
    float **din, **dout;
    sf_file Fin, Fout;

    sf_init(argc,argv);

    Fin = sf_input("in");
    Fout = sf_output("out");

    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");
    if ( !(mode=sf_getstring("mode")) ) mode = "lr";
    /* mode=lr or ud, flip left and right or up and down */

    din = sf_floatalloc2(n1,n2);
    dout = sf_floatalloc2(n1,n2);
    sf_floatread(din[0],n1*n2,Fin);

    if(0==strcmp(mode,"lr")){//flip left and right
    	for(i2=0; i2<n2; i2++)
    	for(i1=0; i1<n1; i1++)
    	{
   		dout[i2][i1]=din[n2-1-i2][i1];
    	}
    }
    if(0==strcmp(mode,"ud")){//flip up and down
    	for(i2=0; i2<n2; i2++)
    	for(i1=0; i1<n1; i1++)
    	{
   		dout[i2][i1]=din[i2][n1-1-i1];
    	}
    }

    sf_floatwrite(dout[0],n1*n2,Fout);

    exit(0);
}

