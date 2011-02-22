/* Take complex comjugate of a (complex) dataset */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
    int n[SF_MAX_DIM], dim, len, i;
    sf_complex *inp, *inpT;
    sf_datatype type;
    sf_file in, out;

    sf_init(argc,argv);
    in   = sf_input("in");
    out  = sf_output("out");

    /* read input dimension */
    dim = sf_filedims(in,n);

    len = 1;
    for (i=0; i < dim; i++) {
	len *= n[i];
    }

    /* read input data type */
    type = sf_gettype(in);

    if (type != SF_COMPLEX) sf_error("input should be complex.");

    inp = sf_complexalloc(len);
    sf_complexread(inp,len,in);

    inpT = sf_complexalloc(len);

    for (i=0; i < len; i++) {
	inpT[i] = sf_cmplx(crealf(inp[i]),-cimagf(inp[i]));
    }

    sf_complexwrite(inpT,len,out);

    exit(0);
}
