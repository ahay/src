/* Creating an edge file for carpet flattening in 2D. */
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
  int n1, i, j, n2, ne, ns, pair[2];
  sf_file inp, edge;

  sf_init(argc, argv);
  inp = sf_input("in");
  edge = sf_output("out");

  if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
  if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");

  if (!sf_getint("shift",&ns)) ns=1;
  /* maximum shift in traces */

  ne = 2*n2*ns - ns*(ns+1);
  
  sf_putint(edge,"n1",2);
  sf_putint(edge,"n2",ne);
  sf_settype (edge, SF_INT);
  
  for (i=0; i < n2; i++) {
    pair[0] = i;
    for (j=0; j < n2; j++) {
      if (i != j && SF_ABS(i-j) <= ns) {
	pair[1] = j;
	sf_intwrite(pair,2,edge);
      }
    }
  }

  exit(0);
}
