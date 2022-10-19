/* Calculate the slope by fitting a line to a set of points in 2-D. */
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

#include <rsf.h>
			     
int main(int argc, char* argv[])
{
    int id, n1, n2, i1, i2;
    float *table=NULL, *trace=NULL;
    float x, y, o1, d1, a, b, sx, sx2, sxy, sy, det;
    sf_file in, out;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    if (!sf_histint(in,"n1",&n1)) sf_error ("Need n1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error ("Need o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error ("Need d1= in input");
    
    n2 = sf_leftsize(in,1);
    
    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);
    
    trace = sf_floatalloc(n1);
    table = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);
	
	sf_floatread(trace,n1,in);
	
	sx = sx2 = sxy = sy = 0.;
	
	for (id=0; id < n1; id++) {
	    x = o1 + id*d1;
	    y = trace[id];
	    sx += x;
	    sx2 += x*x;
	    sxy += x*y;
	    sy += y;
	}
	det = n1*sx2-sx*sx;
	
	if (0.==det) sf_error("zero determinant at trace %d",i2);
	
	a = (n1*sxy-sx*sy)/det;
	b = (sy*sx2-sx*sxy)/det;
	
	sf_warning("b=%f",b);
	
	for (i1=0; i1 < n1; i1++) {
	    x = o1 + i1*d1;
	    table[i1] = a ;
	}

	sf_floatwrite(table,n1,out);
    }
    
    sf_warning(".");
    
    exit(0);
}

