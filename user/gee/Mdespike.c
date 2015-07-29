/* Remove spikes in by sliding 1-D medians. */
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
    int n1, wide, i1, shift, i, i2, n2;
    float *data, *signal, *win;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    data = sf_floatalloc(n1);
    signal = sf_floatalloc(n1);

    if (!sf_getint("wide",&wide)) wide=7;
    /* sliding window width */

    win = sf_floatalloc(wide);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	
	for (i1=0; i1 < n1; i1++) {
	    shift = SF_MAX (0, SF_MIN (n1-wide, i1-wide/2));
	    for (i=0; i < wide; i++) {
		win[i] = data[shift+i];
	    }
	    signal[i1] = sf_quantile(wide/2,wide,win);
	}
	
	sf_floatwrite(signal,n1,out);
    }	

    exit(0);
}

/* 	$Id: Mdespike.c 11860 2014-02-22 16:38:50Z sfomel $	 */
