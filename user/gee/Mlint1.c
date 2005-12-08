/* Linear interpolation */
/*
  Copyright (C) 2005 University of Texas at Austin
  
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

#include "lint1.h"

int main(int argc,char* argv[])
{
    int nm, nd;
    float *m, *d, *c, o1, d1;
    bool adj;
    sf_file mod, dat, coord;

    sf_init(argc,argv);
    coord = sf_input("coord");
    
    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	dat = sf_input("in");
	mod = sf_output("out");
    } else {
	mod = sf_input("in");
	dat = sf_output("out");
    }
    
    if (!sf_histint(coord,"n1",&nd)) sf_error("No n1= in coord");

    if (adj) {
	if (!sf_histint(dat,"n1",&nm) || nm != nd) sf_error("Need n1=%d in input",nd);

	if (!sf_getint("n1",&nm)) sf_error("Need n1=");
	/* regular axis size (for adj=y) */
	if (!sf_getfloat("o1",&o1)) sf_error("Need o1=");
	/* regular axis origin (for adj=y) */
	if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
	/* regular axis sampling (for adj=y) */

	sf_putint(mod,"n1",nm);
	sf_putfloat(mod,"o1",o1);
	sf_putfloat(mod,"d1",d1);
    } else {
	if (!sf_histint(mod,"n1",&nm)) sf_error("No n1= in input");
	if (!sf_histfloat(mod,"o1",&o1)) sf_error("No o1= in input");
	if (!sf_histfloat(mod,"d1",&d1)) sf_error("No d1= in input");

	sf_putint(dat,"n1",nd);
    }

    m = sf_floatalloc(nm);
    d = sf_floatalloc(nd);
    c = sf_floatalloc(nd);

    sf_floatread(c,nd,coord);
    lint1_init(o1,d1,c);

    if (adj) {
	sf_floatread(d,nd,dat);
    } else {
	sf_floatread(m,nm,mod);
    }

    lint1_lop(adj,false,nm,nd,m,d);

    if (adj) {
	sf_floatwrite(m,nm,mod);
    } else {
	sf_floatwrite(d,nd,dat);
    }

    exit(0);
}
