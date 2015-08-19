/* Test k-D tree algorithm. */
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
    int nt, nd;
    float dist;
    float **points, *point;
    kd_node tree, near;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nd)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nt)) sf_error("No n2= in input");
    sf_putint(out,"n2",1);

    points = sf_floatalloc2(nd,nt);
    sf_floatread(points[0],nd*nt,inp);

    tree = kd_tree(points,nt,nd);

    point = sf_floatalloc(nd);
    if (!sf_getfloats("point",point,nd)) sf_error("Need point=");
    
    dist = SF_HUGE;
    kd_nearest(tree, point, 0, nd, &near, &dist);

    sf_floatwrite(kd_coord(near),nd,out);

    exit(0);
}
