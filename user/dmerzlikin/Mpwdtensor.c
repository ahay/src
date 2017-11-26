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

int main (int argc, char *argv[])
{
    int n1,n2,n3, n12,n123, i1,i2,i3, index;
    float *pxpx, *pxpy, *pypy;
    float *l1, *l2, *uv, *uh, *vv, *vh;
    float T, D, norm, swap;
    sf_file in, in2, in3, out, out2, uver, uhor, vver, vhor;
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
    uver = sf_output ("uver"); // file containing vertical component of biggest eigen vector   
    uhor = sf_output ("uhor"); // file containing horizontal component of biggest eigen vector
    vver = sf_output ("vver"); // file containing vertical component of smallest eigen vector   
    vhor = sf_output ("vhor"); // file containing horizontal component of smallest eigen vector
    
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
    uv = sf_floatalloc(n12);
    uh = sf_floatalloc(n12);	
    vv = sf_floatalloc(n12);
    vh = sf_floatalloc(n12);
    
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

			if ( (fabs(pxpy[index])) <= eps) pxpy[index] = 0.0;

			if (pxpy[index] != 0.0){

				/* intermediate variables */
				T = pxpx[index] + pypy[index];
			
				D = pxpx[index]*pypy[index] - pxpy[index]*pxpy[index];

				D = T*T/4.0 - D;

				if (D < 0.0){

					sf_warning("Discriminant is lower than zero - complex roots at [%d,%d,%d] \n",i1,i2,i3);
					D = fabs(D);

				}

				D = sqrt(D);

				/* eigen values */
				l1[index] = 0.5*T + D;
				l2[index] = 0.5*T - D;

				if ( fabs(l2[index]) > fabs(l1[index]) ){

					/* l1 should be the biggest */	
					swap = l1[index];
					l1[index] = l2[index];
					l2[index] = swap;
	
				}

				/* components of the largest eigen value eigen vector */
				uv[index] = l1[index] - pypy[index];
				uh[index] = pxpy[index];

				norm = sqrt(uv[index]*uv[index] + uh[index]*uh[index]);
				
				if (norm != 0.0 && normalize) uv[index] /= norm; uh[index] /= norm;

				/* components of the smallest eigen value eigen vector */
				/*vv[index] = l2[index] - pypy[index];
				vh[index] = pxpy[index];

				norm = sqrt(vv[index]*vv[index] + vh[index]*vh[index]);
				
				if (norm != 0.0 && normalize) vv[index] /= norm; vh[index] /= norm;*/

				if (uv[index] < 0.0){

					uv[index] = -uv[index];
					uh[index] = -uh[index];

				}

				vv[index] = -uh[index];
				vh[index] = uv[index];

				if (vv[index] < 0.0){

					vv[index] = -vv[index];
					vh[index] = -vh[index];

				}
				

			} else {

				/* components of the largest eigen value eigen vector */
				uv[index] = 1.0;
				uh[index] = 0.0;

				/* components of the smallest eigen value eigen vector */
				vv[index] = 0.0;
				vh[index] = 1.0;

				l1[index] = 1.0;
				l2[index] = 1.0;
			
			}		

		}/* done with time loop */ 

	}/* done with iline loop */
	
	/*write all the files*/
	sf_floatwrite(l1,n12,out);
	sf_floatwrite(l2,n12,out2);
	sf_floatwrite(uv,n12,uver);
	sf_floatwrite(uh,n12,uhor);
        sf_floatwrite(vv,n12,vver);
	sf_floatwrite(vh,n12,vhor);

    }/* done with xline loop */
    
    exit (0);
}


