/* 1-D FFT wave extrapolation using Cosine FT */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "lowrank1.h"

int main(int argc, char* argv[]) 
{
    int nx, nt, m1, m2, it, n2;
    float dt, *curr, **lft, **rht, **mid;
    sf_file inp, out, prop, left, right;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1= in input");

    if (!sf_getint("nt",&nt)) sf_error("No nt= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt= in input");

    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.);
    sf_putstring(out,"label2","Time");
    sf_putstring(out,"unit2","s");

    curr = sf_floatalloc(nx);

    prop = sf_input("prop");
    if (SF_FLOAT != sf_gettype(prop)) sf_error("Need float prop");

    left = sf_input("left");
    right = sf_input("right");

    if (!sf_histint(prop,"n1",&m1)) sf_error("No n1= in prop");
    if (!sf_histint(prop,"n2",&m2)) sf_error("No n2= in prop");
    
    if (!sf_histint(left,"n1",&n2) || n2 != nx) sf_error("Need n1=%d in left",nx);
    if (!sf_histint(left,"n2",&n2) || n2 != m1) sf_error("Need n2=%d in left",m1);
	
    if (!sf_histint(right,"n1",&n2) || n2 != m2) sf_error("Need n1=%d in right",m2);
    if (!sf_histint(right,"n2",&n2) || n2 != nx) sf_error("Need n2=%d in right",nx);
    
    lft = sf_floatalloc2(nx,m1);
    mid = sf_floatalloc2(m1,m2);
    rht = sf_floatalloc2(m2,nx);

    sf_floatread(lft[0],nx*m1,left);
    sf_floatread(rht[0],nx*m2,right);
    sf_floatread(mid[0],m1*m2,prop);

    sf_fileclose(left);
    sf_fileclose(right);
    sf_fileclose(prop);

    lowrank1_init(nx,m1);
    sf_floatread (curr,nx,inp);

    /* propagation in time */
    for (it=0; it < nt; it++) {
	lowrank1_step(m1,m2,lft,mid,rht,curr);
	sf_floatwrite(curr,nx,out);
    }
    
    exit(0);
}
