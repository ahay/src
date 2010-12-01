/* Diffraction focusing test. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
    int nt, nx;
    float v0, v1;
    float **data;
    sf_file inp, out;

    inp = sf_input("in");
    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_getfloat("v0",&v0)) v0=SF_EPS;
    /* initial velocity */

    if (!sf_getfloat("v",&v1)) sf_error("Need v=");
    /* final velocity */

    data = sf_floatalloc2(nt,nx);

    sf_floatread(data[0],nt*nx,inp);

    sf_floatwrite(data[0],nt*nx,out);

    exit(0);
}
