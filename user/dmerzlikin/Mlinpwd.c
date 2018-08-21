/* Chain of 3D Path Integral, Azimuthal Plane-Wave Destruction and Kirchhoff migration (based on sfmig3)

works only for zero offset

make sure nh = 1 dh = 1.0 h0 = 0.0 offset file is not used
 
there are flags to disable Azimuthal PWD (Plane-Wave Destruction), P (Path-Integral Filter) and L (Kirchhoff modelling/migration)

no regularization

can be expressed for forward as:

data = P PWD L ( diffractions ) or as a matrix

| P PWD L   P PWD L | | diffractions | = | data |                     

can be expressed for adjoint as:

adjoint diffractions = L^T PWD^T P^T data or as a matrix

| diffractions | = | L^T PWD^T P^T | | data |

*/
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
#include <math.h>
#include "allp3.h"

int main(int argc, char* argv[])
{
    bool sm, adj;
    int nt,nx,ny,n12,n123, nw,nj1,nj2, n4;
    float dt,dx,dy, ot,ox,oy;
    sf_file in, dip, out;
    float *pp1, *pp2, *data, *pwddata;
    allpass ap, aq;

    sf_init(argc,argv);
    in = sf_input ("in");
    dip = sf_input ("dip");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(dip)) sf_error("Need float type");
   
    /* get dimensions from input */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in in");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in in");
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in in");
    if (!sf_histint(in,"n4",&n4)) n4=1;
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in in");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in in");
    if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in in");
    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in in");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in in");
    if (!sf_histfloat(in,"o3",&oy)) sf_error("No o3= in in");

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
	
    if (nw < 1 || nw > 3) sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing iline */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing xline */

    if (!sf_getbool("sm",&sm)) sm=true;
    /* if perform AzPWD filtering */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    n12 = nt*nx;
    n123 = n12*ny;

    if (adj == false){/* if fwd - two volumes: PWDx and PWDy */

    	sf_putint(out,"n4",2);
    	n4=2;

    } else {/* if adj - one volume */

    	sf_putint(out,"n4",1);
    	n4=1;

    }

    /* allocate space for data */
    data = sf_floatalloc(n123); 
    pwddata = sf_floatalloc(2*n123); 

    /* allocate space for dip */
    pp1 = sf_floatalloc(n123);
    pp2 = sf_floatalloc(n123);

    /* reading */
    if (adj == false){
    	sf_floatread(data,n123,in);
    } else {
	sf_floatread(pwddata,2*n123,in);
    }
    /* iline dip */
    sf_floatread(pp1,n123,dip);
    /* xline dip */
    sf_floatread(pp2,n123,dip);

    /* initialize linear PWD filter */

    /* iline */
    ap = allpass_init(nw, nj1, nt,nx,ny, pp1);

    /* xline */
    aq = allpass_init(nw, nj2, nt,nx,ny, pp2);

    /* iline + xline */
    allpass3_init(ap,aq);

    if (sm){/* perform AzPWD */

	allpass3_lop(adj,false,n123,2*n123,data,pwddata);

    }/* AzPWD */

    if (adj == false){
    	sf_floatwrite(pwddata,n4*n123,out);
    } else {
        sf_floatwrite(data,n4*n123,out);
    }

    exit(0);

}
