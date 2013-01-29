/* Automatic event PICKing */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>


int main(int argc, char* argv[])
{
    int n1, n2, n3, i1, i2 ,i3;
    float o1, d1;
    float **u1, *u2, c=0.0;

    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");  	/* seismic age */
    out = sf_output ("out"); 	/*  */


    if (!sf_histint(in, "n1", &n1)) sf_error("n1 error");
    if (!sf_histint(in, "n2", &n2)) sf_error("n2 error");
    if (!sf_histint(in, "n3", &n3)) n3=1;
    if (!sf_histfloat(in, "o1", &o1)) sf_error("o1");
    if (!sf_histfloat(in, "d1", &d1)) sf_error("d1");
    sf_unshiftdim(in,out,1);

    u1 = sf_floatalloc2(n1, n2);
    u2 = sf_floatalloc(n2);

    for(i3=0; i3<n3; i3++)
    {
	sf_floatread(u1[0], n1*n2, in);
	for(i2=0; i2<n2; i2++)
	{
	    for(i1=0; i1<n1; i1++)
		if(u1[i2][i1] >= c) break;
	    u2[i2] = d1*i1+o1;
	}
	sf_floatwrite(u2, n2, out);
    }

    free(u1[0]);
    free(u1);
    free(u2);

    return 0;
}



