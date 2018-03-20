/* structure tensor estimation based on plane wave destruction. */

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

#include "eigen1.h"

int main (int argc, char *argv[])
{
    int n1,n2,n3, n12,n123, i1,i2,i3, index;
    float *pxpx, *pxpy, *pypy;
    float *l1, *l2, *ux, *uy, *vx, *vy;
    float T, D, norm, swap;
    sf_file in, in2, in3, out, out2, uxf, uyf, vxf, vyf;
    float eps; // eps tolerance
    bool normalize;// normalize eigen vectors
    // add flags so that it outputs only the perpendicular to structures component
    // otherwise you will have lots of files all the time

    sf_init(argc,argv);
    in = sf_input ("in"); // <pwdx * pwdx> - where pwdx/y - pwd applied in inline/crossline direction
    in2 = sf_input ("in2"); // <pwdx * pwdy> - and <> denotes structure oriented smoothing
    in3 = sf_input ("in3"); // <pwdy * pwdy> - with the same dip distribution
    out = sf_output ("out"); // file containing eigen values of structure tensor (only the biggest one)
    out2 = sf_output ("out2"); // file containing eigen values of structure tensor (only the smallest one)
    uxf = sf_output ("ux"); // file first component of biggest eigen vector   
    uyf = sf_output ("uy"); // file second component of biggest eigen vector
    vxf = sf_output ("vx"); // file first component of smallest eigen vector   
    vyf = sf_output ("vy"); // file second component of smallest eigen vector
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    /* reading files dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histint(in,"n3",&n3)) n3=1;
    n12 = n1*n2;
    n123 = n12*n3;

    if (!sf_getfloat("eps",&eps)) eps=0.00001; // tolerance
    if (!sf_getbool("normalize",&normalize)) normalize=false; // normalize eigen vectors

    /* memory allocation */
    pxpx = sf_floatalloc(n12);
    pxpy = sf_floatalloc(n12);
    pypy = sf_floatalloc(n12);
    
    l1 = sf_floatalloc(n12);
    l2 = sf_floatalloc(n12);
    ux = sf_floatalloc(n12);
    uy = sf_floatalloc(n12);	
    vx = sf_floatalloc(n12);
    vy = sf_floatalloc(n12);
    
    /* xline loop */
    for (i3=0; i3<n3; i3++){
	
	/*read all the files*/
	sf_floatread(pxpx,n12,in);
        sf_floatread(pxpy,n12,in2);
        sf_floatread(pypy,n12,in3);

	/* iline loop */
	for (i2=0; i2<n2; i2++){
		
		/* time  loop */
		for (i1=0; i1<n1; i1++){

			index = i2*n1 + i1;

			//sf_warning("be4 ux=%f uy=%f vx=%f vy=%f",ux[index],uy[index],vx[index],vy[index]);

			solveSymmetric22(pxpx[index], pxpy[index], pypy[index], &ux[index],&uy[index],&vx[index],&vy[index],&l1[index],&l2[index]);

			if (ux[index] < 0.0f){

				ux[index] = -ux[index];

				uy[index] = -uy[index];

			}

			vx[index] = -uy[index];

			vy[index] = ux[index];

			if(l2[index]<0.0f) l2[index]=0.0f;

			if(l1[index]<l2[index]) l1[index]=l2[index];

			//sf_warning("assigned ux=%f uy=%f vx=%f vy=%f",ux[index],uy[index],vx[index],vy[index]);

		}/* done with time loop */ 

	}/* done with iline loop */
	
	/*write all the files*/
	sf_floatwrite(l1,n12,out);
	sf_floatwrite(l2,n12,out2);
	sf_floatwrite(ux,n12,uxf);
	sf_floatwrite(uy,n12,uyf);
        sf_floatwrite(vx,n12,vxf);
	sf_floatwrite(vy,n12,vyf);

    }/* done with xline loop */
    
    exit (0);
}


