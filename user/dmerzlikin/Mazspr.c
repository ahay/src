/* Combining Sprays: Simply Input Sprays in In-line And Cross-line */
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
    bool sm;
    int nt,nx,ny,n12,n123, nw,nj1,nj2, n4, i;
    float dt,dx,dy, ot,ox,oy;
    float conv;
    sf_file in, in2, azin, out;
    float *data, *data2, *az, *comb;

    sf_init(argc,argv);
    in = sf_input ("in");/* 1 - spray_x */
    in2 = sf_input ("spry");/* 1 - spray_y */
    azin = sf_input ("az");/* azimuth */
    out = sf_output ("out");
   
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

    //if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
	
    //if (nw < 1 || nw > 3) sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    //if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing iline */
    //if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing xline */

    if (!sf_getbool("sm",&sm)) sm=true;
    /* if perform AzPWD filtering */

    //if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    n12 = nt*nx;
    n123 = n12*ny;

    /* allocate space for data */
    data = sf_floatalloc(n123);
    data2 = sf_floatalloc(n123);
    comb = sf_floatalloc(n123);

    /* allocate space for azimuth */
    az = sf_floatalloc(n123);

    /* reading data: spray_x */
    sf_floatread(data,n123,in);

    /* reading data: spray_x */
    sf_floatread(data2,n123,in2);

    /* reading azimuth */
    sf_floatread(az,n123,azin);

    conv = 3.14159265 / 180;

    if (sm) {

	for (i=0; i<n123; i++){	

		data[i] = data[i]*(cosf(conv*az[i]));

		data2[i] = data2[i]*(sinf(conv*az[i]));			

		comb[i] = data[i] + data2[i];

	}

    }

    sf_floatwrite(comb,n123,out);

    exit(0);

}
