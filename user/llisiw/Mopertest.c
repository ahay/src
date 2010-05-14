/* First-arrival tomography. */
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

#include "upgrad.h"

int main(int argc, char* argv[])
{
    bool inv;
    int dim=2, n[2];
    float o[2], d[2], *t, *x, *r;
    upgrad upg;
    sf_file in, out, time;

    sf_init(argc,argv);
    in   = sf_input("in");
    time = sf_input("time");
    out  = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;

    if (!sf_histint(in,"n1",&n[0]))   sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n[1]))   sf_error("No n2= in input");
    if (!sf_histfloat(in,"o1",&o[0])) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d[1])) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o2",&o[0])) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&d[1])) sf_error("No d2= in input");

    t = sf_floatalloc(n[0]*n[1]);
    x = sf_floatalloc(n[0]*n[1]);
    r = sf_floatalloc(n[0]*n[1]);
    sf_floatread(t,n[0]*n[1],time);
    sf_floatread(x,n[0]*n[1],in);

    upg = upgrad_init(dim,n,d);

    upgrad_set(upg,t);

    if (inv) {
	upgrad_solve(upg,x,r,NULL);
    } else {
	upgrad_forw(upg,x,r);
    }

    sf_floatwrite(r,n[0]*n[1],out);

    exit(0);
}
