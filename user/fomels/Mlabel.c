/* Connected-component labeling */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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

#include "label.h"

int main (int argc, char* argv[])
{
    int i1, i2, n1, n2, maxl, a,c,d, **img, **lbl;
    sf_file inp, out;

    sf_init (argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_INT != sf_gettype(inp)) sf_error("Need int type in input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");

    img = sf_intalloc2(n1,n2);
    lbl = sf_intalloc2(n1,n2);
    
    sf_intread(img[0],n1*n2,inp);

    maxl = 0;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (img[i2][i1]) maxl++;
	    lbl[i2][i1]=0;
	}
    }
    
    label_init(maxl);

    	    
    /* -------------
       | a | b | c |
       -------------
       | d | e |   |
       -------------
       |   |   |   |
       ------------- */
    
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (!img[i2][i1]) continue; /* background */
	    
	    if (i1 > 0 && img[i2][i1-1]) {
		lbl[i2][i1] = lbl[i2][i1-1];
	    } else if (i2+1 < n2 && i1 > 0 && img[i2+1][i1-1]) {
		c = lbl[i2+1][i1-1];
		lbl[i2][i1] = c;
		
		if (i2 > 0 && img[i2-1][i1-1]) {
		    a = lbl[i2-1][i1-1];
		    label_union(c,a);
		}
	    
		if (i2 >0 && img[i2-1][i1]) {
		    d = lbl[i2-1][i1];
		    label_union(c,d);
		}
	    } else if (i2 > 0 && i1 > 0 && img[i2-1][i1-1]) {
		lbl[i2][i1] = lbl[i2-1][i1-1];
	    } else if (i2 >0 && img[i2-1][i1]) {
		lbl[i2][i1] = lbl[i2-1][i1];
	    } else {
		lbl[i2][i1] = label_new();
	    }
	}
    }

    label_flatten();

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    img[i2][i1] = label_find(lbl[i2][i1]);
	}
    }

    sf_intwrite(img[0],n1*n2,out);
		        
    exit(0);
}
