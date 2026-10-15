/* Creating an edge file for carpet flattening in 3D. */
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
  int n1, n2, n3, i2, i3, j2, j3, ne, ns, pair[2];
  sf_file inp, edge;

  sf_init(argc, argv);
  inp = sf_input("in");
  edge = sf_output("out");

  if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
  if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
  if (!sf_histint(inp,"n3",&n3)) n3=1;

  if (!sf_getint("shift",&ns)) ns=1;
  /* maximum shift in traces */

  ne = 0;
  for (i3=0; i3 < n3; i3++) {
    for (i2=0; i2 < n2; i2++) {
      for (j3=SF_MAX(0,i3-ns); j3 <= SF_MIN(n3-1,i3+ns); j3++) {
	for (j2=SF_MAX(0,i2-ns); j2 <= SF_MIN(n2-1,i2+ns); j2++) {
	  if (i2 != j2 || i3 != j3) ne++;
	}
      }
    }
  }

  sf_putint(edge,"n1",2);
  sf_putint(edge,"n2",ne);
  sf_putint(edge,"n3",1);
  sf_settype (edge, SF_INT);


  for (i3=0; i3 < n3; i3++) {
    for (i2=0; i2 < n2; i2++) {
      pair[0] = i2+n2*i3;
      for (j3=SF_MAX(0,i3-ns); j3 <= SF_MIN(n3-1,i3+ns); j3++) {
	for (j2=SF_MAX(0,i2-ns); j2 <= SF_MIN(n2-1,i2+ns); j2++) {
	  if (i2 != j2 || i3 != j3) {
	    pair[1] = j2+n2*j3;
	    sf_intwrite(pair,2,edge);
	  }
	}
      }
    }
  }
  
  exit(0);
}
