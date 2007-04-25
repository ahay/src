/* Remove spikes in by sliding 3-D medians. */
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

int main(int argc, char* argv[]) 
{
    int wide, shift1, shift2, shift3, i, j, k, i1, n1, i2, n2, i3, n3, i4, n4;
    float ***data, ***signal, ***win;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n4 = sf_leftsize(in,3);

    data = sf_floatalloc3(n1,n2,n3);
    signal = sf_floatalloc3(n1,n2,n3);

    if (!sf_getint("wide",&wide)) wide=5;
    /* sliding window width */

    win = sf_floatalloc3(wide,wide,wide);

    for (i4=0; i4 < n4; i4++) {
	sf_floatread(data[0][0],n1*n2*n3,in);	
	for (i3=0; i3 < n3; i3++) { shift3 = SF_MAX (0, SF_MIN (n3-wide, i3-wide/2 - 1));
	for (i2=0; i2 < n2; i2++) { shift2 = SF_MAX (0, SF_MIN (n2-wide, i2-wide/2 - 1));
	for (i1=0; i1 < n1; i1++) { shift1 = SF_MAX (0, SF_MIN (n1-wide, i1-wide/2 - 1));
	for (i=0; i < wide; i++) {
	    for (j=0; j < wide; j++) {
		for (k=0; k < wide; k++) {
		    win[i][j][k] = data[shift3+k][shift2+i][shift1+j];
		}
	    }
	}
	signal[i3][i2][i1] = sf_quantile(wide*wide*wide/2,wide*wide*wide,win[0][0]);
	}}}
	
	sf_floatwrite(signal[0][0],n1*n2*n3,out);
    }	

    exit(0);
}

/* 	$Id$	 */
