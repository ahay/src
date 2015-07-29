/* First-arrival Traveltime Tomography (linear operator) */
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

#include "ftoper.h"

int main(int argc, char* argv[])
{
    bool adj;
    int n[3], nt, ns, **m, is, k;
    float d[3], **t0, *dt, *ds;
    sf_file in, out, time, mask;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */    

    if (NULL == sf_getstring("time")) sf_error("Need time=");

    time = sf_input("time");

    if (!sf_histint(time,"n1",&n[0])) sf_error("No n1= in input");
    if (!sf_histint(time,"n2",&n[1])) sf_error("No n2= in input");
    if (!sf_histint(time,"n3",&n[2])) n[2]=n[1];

    if (!sf_histfloat(time,"d1",&d[0])) sf_error("No n1= in input");
    if (!sf_histfloat(time,"d2",&d[1])) sf_error("No n2= in input");
    if (!sf_histfloat(time,"d3",&d[2])) d[2]=d[1];

    if (!sf_histint(time,"n4",&ns)) ns=1;

    nt = n[0]*n[1]*n[2];

    t0 = sf_floatalloc2(nt,ns);
    sf_floatread(t0[0],nt*ns,time);

    m = sf_intalloc2(nt/n[0],ns);

    if (NULL == sf_getstring("mask")) {
	mask = NULL;

	for (is=0; is < ns; is++)
	    for (k=0; k < n[1]*n[2]; k++)
		m[is][k] = 1;
    } else {
	mask = sf_input("mask");
	sf_intread(m[0],nt/n[0],mask);
	sf_fileclose(mask);
    }
    
    ftoper_init(n,d,ns,m);

    ftoper_set(t0);

    dt = sf_floatalloc(ns*nt);
    ds = sf_floatalloc(nt);

    if (adj) {
	sf_floatread(dt,ns*nt,in);
	sf_putint(out,"n4",1);
    } else {
	sf_floatread(ds,nt,in);	
	sf_putint(out,"n4",ns);
    }
    
    if (adj) {
	ftoper_oper(true,false,nt,ns*nt,ds,dt);
    } else {
	ftoper_oper(false,false,nt,ns*nt,ds,dt);	    
    }

    if (adj) {
	sf_floatwrite(ds,nt,out);
    } else {
	sf_floatwrite(dt,ns*nt,out);
    }

    exit(0);
}
