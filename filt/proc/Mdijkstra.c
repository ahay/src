/* Dijkstra shortest-path algorithm in 2-D */
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

#include "dijkstra.h"

int main (int argc, char *argv[])
{
    int n1,n2,ref1,ref2,n12;
    float **p,**q;
    sf_file out, cost;

    sf_init(argc,argv);
    cost = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(cost,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(cost,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    if (!sf_getint("ref1",&ref1)) ref1=0;
    if (!sf_getint("ref2",&ref2)) ref2=0;
    /* reference trace */

    p = sf_floatalloc2(n1,n2);
    q = sf_floatalloc2(n1,n2);

    sf_floatread(p[0],n12,cost);
    sf_floatread(q[0],n12,cost);

    sf_putint(out,"n3",1);

    dijkstra_init(n1,n2);
    dijkstra(ref1,ref2,p,q);

    sf_floatwrite(dijsktra_cost(),n12,out);

    exit(0);
}
