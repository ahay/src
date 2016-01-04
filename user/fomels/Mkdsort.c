/* Sort entries based on k-D tree. */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "kdtree.h"

int main(int argc, char* argv[])
{
    int n, nt, nd, i3, n3, *order;
    float **points;
    kd_node tree;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nd)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input");
    n3 = sf_unshiftdim(inp,out,1);
    sf_putint(out,"n1",nt);
    sf_settype(out,SF_INT);

    points = sf_floatalloc2(nd,nt);
    order = sf_intalloc(nt);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(points[0],nd*nt,inp);

	tree = kd_tree(points,nt,nd);
	n = kd_order(tree,order,0);
	if (nt != n) sf_error("Need %d, got %d",nt,n);

	sf_intwrite(order,nt,out);

	free_tree(tree,nt);
    }

    exit(0);
}
