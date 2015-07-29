/* Post-stack velocity continuation by implicit finite differences */
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

#include "velcon.h"

int main(int argc, char* argv[])
{
    int n1, n2, n12, nv, inv;
    bool adj, add;
    float *modl, *data, o1,d1,o2,d2, vel,v0;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

    if (!sf_getfloat("vel",&vel)) vel=0.75;
    /* final velocity */
    if (!sf_getfloat("v0",&v0)) v0=0.;
    /* starting velocity */
    if (!sf_getint("nv",&nv)) nv=n1;
    /* number of steps */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    if (!sf_getbool("add",&add)) add=false;
    /* addition flag */

    if (!sf_getint("inv",&inv)) inv=0;
    /* amplitude type */

    modl = sf_floatalloc(n12);
    data = sf_floatalloc(n12);

    sf_floatread (adj? data: modl,n12,in);

    velcon_init (inv,vel,v0,o1,n1,n2,nv,d1,d2);
    velcon_lop (adj, add, n12, n12, modl, data);

    sf_floatwrite (adj? modl: data,n12,out);

    exit(0);
}
