/* Stack with antialiasing */
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
    int nt, nx, nv, ix, it;
    bool adj, inv, box;
    float vel, antialias, o1, o2, d1, d2, slow, t;
    float **gather, *stack, *coord, *delt, *amp;
    sf_file inp, out;

    sf_init (argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (!sf_getbool("inv",&inv)) inv=false;
    /* inverse flag */

    if (adj==inv) {
	if (!sf_histint(inp,"n2",&nv)) sf_error("No n2= in input");
	if (!sf_getint("n2",&nx)) nx=15;
	sf_putint(out,"n2",nx);
    } else {
	if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_getint("n2",&nv)) nv=1;
	sf_putint(out,"n2",nv);
    }

    if (!sf_histfloat(inp,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(inp,"d2",&d2)) sf_error("No d2= in input");

    if (!sf_getfloat("vel",&vel)) vel=1.5;
    /* velocity */

    if (!sf_getfloat("antialias",&antialias)) antialias=1.;
    /* antialiasing */

    if (!sf_getbool("box",&box)) box=false;
    /* box antialiasing */

    gather = sf_floatalloc2(nt,nx);
    stack = sf_floatalloc(nt);
    coord = sf_floatalloc(nt);
    delt = sf_floatalloc(nt);
    amp = sf_floatalloc(nt);

    if (adj) {
	sf_floatread(gather[0],nt*nx,inp); 
    } else {
	sf_floatread (stack,nt,inp); 
    }

    sf_adjnull(adj,false,nt,nt*nx,stack,gather[0]);

    sf_aastretch_init (box, nt, o1, d1, nt);

    for (ix=0; ix < nx; ix++) {
	slow = (o2 + d2*ix)/vel;

	for (it=0; it < nt; it++) {
	    t = o1 + it*d1;
	    coord[it] = hypotf(t,slow);
	    if (coord[it] > 0.) {
		delt[it] = fabsf(antialias*slow*d2/(vel*coord[it]));
		amp[it] = (t / coord[it]) / sqrtf (coord[it]);
	    } else {
		delt[it] = 0.;
		amp[it] = 0.;
	    }
	}

	sf_aastretch_define (coord, delt, NULL);
	sf_aastretch_lop (adj,true,nt,nt,stack,gather[ix]);
    }

    if (adj) {
	sf_floatwrite(stack,nt,out);
    } else {
	sf_floatwrite(gather[0],nt*nx,out);
    }

    exit(0);
}

