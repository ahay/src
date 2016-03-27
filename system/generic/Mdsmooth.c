/* Multi-dimensional triangle smoothing - derivative. */
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

int main (int argc, char* argv[]) 
{
    int dim, dim1, i, ider, nrep, i2, n1, n2, nderiv;
    int n[SF_MAX_DIM], rect[SF_MAX_DIM];
    char key[6];
    float* data;
    sf_file in, out;

    sf_init (argc, argv);
    in  = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) dim1 = i;
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

    data = sf_floatalloc (n1);

    if (!sf_getint("ider",&ider)) ider=0;
    /* direction of the derivative (0 means no derivative) */
    
    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat smoothing several times */

    if (!sf_getint("nderiv",&nderiv)) nderiv=6;
    /* derivative filter accuracy */

    sf_dtrianglen_init(dim1+1,rect,n);
    
    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	sf_dtrianglen(ider, nrep, nderiv, data);
	
	sf_floatwrite(data,n1,out);
    }    

    exit (0);
}

