/* Convert RSF to velocity picks */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

int main (int argc, char* argv[])
{
    int n1,n2,n3 = 0, i1,i2,i3,skip;
    float **xy, o1,d1,o2,d2;
    sf_file in;

    sf_init (argc,argv);
    in = sf_input("in");
    //out = sf_output("out");

    /* read data */

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");
    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("Need n2=");
    if (!sf_histint(in,"n3",&n2)) n3=1;
   
    if (!sf_getint ("skip",&skip)) skip=4;
    /* number of distance bins */

    if (!sf_histfloat(in,"o1",&o1)) sf_error("Need o1=");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("Need o2=");
//    if (!sf_histfloat(in,"o3",&n1)) o3=o2;

    if (!sf_histfloat(in,"d1",&d1)) sf_error("Need d1=");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("Need d2=");
//    if (!sf_histfloat(in,"d3",&n1)) d3=d2;

    xy = sf_floatalloc2(n1,n2);
//    sf_floatread(xy[0][0],n1*n2*n3,in);

    /* loop through data pairs */
    for (i3=0; i3 < n3; i3++) {
	sf_floatread(xy[0],n1*n2,in);	
	for (i2=0; i2 < n2; i2++) {
            printf( "%-8s%-8d\n", "VFUNC", i2);
	    for (i1=0; i1 < n1; i1+=skip) {
		if (i1%(4*skip)==0 && i1!=0) printf ("\n");
                printf ("%8d%8.0f", i1*1000,xy[i2][i1]);
            }
            printf ("\n");
        }
    }
//    sf_floatwrite (vari,nx,out);
    exit(0);
}

