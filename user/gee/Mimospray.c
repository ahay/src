/* Inversion of constant-velocity nearest-neighbor inverse NMO. */
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
#include "imospray.h"

int main(int argc, char* argv[])
{
    int n1,n2,n12;
    bool adj;
    float d1,d2,o1,o2,v, *model=NULL, *dat=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (adj) {
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");

	sf_putint(out,"n2",1);
    } else {
	if (!sf_getint("n2",&n2)) n2=20;
	/* number of offsets (if inv=n) */
	if (!sf_getfloat("d2",&d2)) d2=200.;
	/* offset sampling (if inv=n) */
	if (!sf_getfloat("o2",&o2)) o2=0.;
	/* offset origin (if inv=n) */

	sf_putint(out,"n2",n2);
	sf_putfloat(out,"d2",d2);
	sf_putfloat(out,"o2",o2);
    }
    n12 = n1*n2;

    if (!sf_getfloat ("v",&v)) v=1000.;
    /* velocity */

    dat = sf_floatalloc(n12);
    model =  sf_floatalloc(n1);

    imospray_init (1./v, o2,d2, o1,d1, n1,n2);

    if (adj) {
	sf_floatread(dat,n12,in);
	imospray_lop (true,false,n1,n12,model,dat);
	sf_floatwrite(model,n1,out);
    } else {	
	sf_floatread(model,n1,in);
	imospray_lop (false,false,n1,n12,model,dat);
	sf_floatwrite(dat,n12,out);
    }

    exit(0);
}
