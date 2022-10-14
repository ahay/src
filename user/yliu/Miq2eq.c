/* Convert interval Q value to equivalent Q value */
/*
  Copyright (C) 2022 Jilin University
  
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

int main (int argc,char *argv[])
{
    bool verb;

    float o1, d1;
    int n1, n2, i1, i2, i;    
    float detf, tp;
    float *teq,*tq;
    sf_file inp, out;

    sf_init (argc,argv);

    inp = sf_input("in");
    out = sf_output("out");    

    if (!sf_getbool("verb",&verb)) verb=false;
	
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&o1)) sf_error("No o1= in input");

    n2=sf_leftsize(inp,1);

    tq = sf_floatalloc(n1);    
    teq = sf_floatalloc(n1);

    
    for (i2=0; i2<n2; i2++) {
	sf_floatread(tq,n1,inp);
	
	for (i1=2; i1<n1; i1++) {
	    tp=0;
	    detf = 0.;
	    tp=i1;
	    for (i=1; i<=tp; i++) {
		
		detf += d1/tq[i] ;
	    }	    
	    teq[i1] = i1*d1/detf;
	}

	for (i1=0; i1<=1; i1++) {
	    teq[i1] = tq[i1];
	}
	
	sf_floatwrite(teq,n1,out);
    }
}
