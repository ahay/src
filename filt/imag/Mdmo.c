/* Kirchhoff DMO with antialiasing by reparameterization. */
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

#include "dmo.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, n12, i3,mint,n,type;
    float *dat1, *dat2, t0,dt,dx,h,velhalf;
    bool adj, inv;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint ("mint",&mint)) mint=2;
    /* starting time sample */

    if (!sf_getint ("n",&n)) n=32;
    /* number of offset samples */

    if (!sf_getbool ("adj",&adj)) adj=false;
    /* adjoint flag */

    if (!sf_getbool ("inv",&inv)) inv=false;
    /* inversion flag */

    if (!sf_getint ("type",&type)) type=1;
    /* type of amplitude (0,1,2,3) */

    if (!sf_getfloat ("h",&h)) sf_error("Need h=");
    /* half-offset */

    if (!sf_getfloat ("velhalf",&velhalf)) velhalf=0.75;
    /* half-velocity */

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

    dat1 = sf_floatalloc(n12);
    dat2 = sf_floatalloc(n12);

    dmo_init (velhalf,inv,t0,dt,dx,n1,n2,h,mint,n,type);

    for (i3=0; i3 < n3; i3++) { 
	sf_floatread(adj? dat2: dat1,n12,in); 

	dmo_lop(adj,false,n12,n12,dat1,dat2);

        sf_floatwrite(adj? dat1: dat2,n12,out);
    }

    exit(0);
}

