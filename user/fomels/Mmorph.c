/* Morphological operations on binary images */
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

int main (int argc, char* argv[])
{
    int i1, i2, j1, j2, n1, n2, **img1, **img2, a;
    char *what;
    sf_file inp, out;
    
    sf_init (argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_INT != sf_gettype(inp)) sf_error("Need uchar type in input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");

    if (NULL == (what = sf_getstring("what"))) what="none";

    img1 = sf_intalloc2(n1,n2);
    img2 = sf_intalloc2(n1,n2);

    sf_intread(img1[0],n1*n2,inp);
	   
    for(i2=1; i2 < n2-1; i2++){
	for(i1=1; i1 < n1-1; i1++){
	    switch (what[0]) {
		case 'd': /* dilation */
		case 'c': /* closing */
		    a = 0;
		    for(j2=i2-1; j2<=i2+1; j2++){
			for(j1=i1-1; j1<=i1+1; j1++){
			    if(img1[j2][j1] > a)
				a = img1[j2][j1];
			}  
		    }  
		    break;
		case 'e': /* erosion */
		case 'o': /* opening */
		    a = 1;
		    for(j2=i2-1; j2<=i2+1; j2++){
			for(j1=i1-1; j1<=i1+1; j1++){
			    if(img1[j2][j1] < a)
				a = img1[j2][j1];
			}  
		    }  
		    break;
		case 'n': /* none */
		default:
		    a = img1[i2][i1];
		    break;
	    }
	    img2[i2][i1] = a;
	}  /* ends loop over j */
    }  /* ends loop over i */

    for(i2=1; i2 < n2-1; i2++){
	for(i1=1; i1 < n1-1; i1++){
	    switch (what[0]) {
		case 'o': /* opening */
		    a = 0;
		    for(j2=i2-1; j2<=i2+1; j2++){
			for(j1=i1-1; j1<=i1+1; j1++){
			    if(img2[j2][j1] > a)
				a = img2[j2][j1];
			}  
		    }  
		    break;
		case 'c': /* closing */
		    a = 1;
		    for(j2=i2-1; j2<=i2+1; j2++){
			for(j1=i1-1; j1<=i1+1; j1++){
			    if(img2[j2][j1] < a)
				a = img2[j2][j1];
			}  
		    }  
		    break;
		default:
		    a = img2[i2][i1];
		    break;
	    }
	    img1[i2][i1] = a;
	}
    }

    sf_intwrite(img1[0],n1*n2,out);

    exit(0);
}
