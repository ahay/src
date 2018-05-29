/* Azimuthal Plane-Wave Destruction */
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
#include "azpwd.h"

int main(int argc, char* argv[])
{
    bool sm, adj;
    int nt,nx,ny,n12,n123, nw,nj1,nj2, n4;
    float dt,dx,dy, ot,ox,oy;
    sf_file in, dip, azin, out;
    float *pp1, *pp2, *data, *pwddata, *az;

    sf_init(argc,argv);
    in = sf_input ("in");
    dip = sf_input ("dip");
    azin = sf_input ("az");
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

    /* allocate space for data */
    data = sf_floatalloc(n123);
    pwddata = sf_floatalloc(n123); 

    /* allocate space for dip */
    pp1 = sf_floatalloc(n123);
    pp2 = sf_floatalloc(n123);

    /* allocate space for azimuth */
    az = sf_floatalloc(n123);

    /* reading data */
    sf_floatread(data,n123,in);
    
    /* reading iline dip */
    sf_floatread(pp1,n123,dip);
    
    /* reading xline dip */
    sf_floatread(pp2,n123,dip);

    /* reading azimuth */
    sf_floatread(az,n123,azin);

    azpwd_init(nt, nx, ny, dt, dx, dy, ot, ox, oy, 
		 nw /* [1,2,3] accuracy order */,
		 pp1, pp2 /* dip distributions */,
		 az /* azimuth distribution */,
                 nj1, nj2 /* antialiasing */);

    if (sm) {

	if (!adj){

		azpwd_lop (adj,false,n123,n123,data,pwddata);

	} else {

		azpwd_lop (adj,false,n123,n123,pwddata,data);

	}   

    }

    sf_floatwrite(pwddata,n123,out);

    exit(0);

}
