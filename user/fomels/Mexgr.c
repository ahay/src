/* Exact group velocity in VTI media */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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
    const int na=361;
    int ia;
    float a, b, a2, da, theta, vel;
    float A, C, B, D;
    sf_complex *plot;
    sf_file out;

    sf_init (argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_complex");

    sf_putint(out,"n1",na);

    if (!sf_getfloat("aq",&A)) A=14.47;
    if (!sf_getfloat("bq",&B)) B=2.28;
    if (!sf_getfloat("cq",&C)) C=9.57;
    if (!sf_getfloat("dq",&D)) D=4.51;

    plot = sf_complexalloc(na);
    da = SF_PI/180.;

    for (ia=0; ia < na; ia++) {
	a = cosf(da*ia);
	b = sinf(da*ia);

	a2 = a*a;
	
	theta = (a2*powf(-(A*(-1 + a2)*(B - C)) + a2*powf(B - C,2) +
			 (1 - a2)*(B*(B + C) + 4*B*D + 2*powf(D,2)) +
			 (B + C)*sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
				       4*(-1 + a2)*a2*powf(B + D,2)),2))/
	    (2.*(-(powf(A,4)*powf(-1 + a2,3)) +
		 powf(a,6)*powf(B - C,2)*(powf(B,2) + powf(C,2)) + powf(-1 + a2,2)*powf(B,3)*
		 (B - a2*B + sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) - 4*(-1 + a2)*a2*powf(B + D,2))) +
		 powf(A,3)*powf(-1 + a2,2)*(-2*B + 4*a2*B - 2*a2*C + 
					    sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
						  4*(-1 + a2)*a2*powf(B + D,2))) +
		 powf(A,2)*(1 - a2)*(powf(a,4)*powf(B,2) +
				     2*powf(-1 + a2,2)*powf(B,2) - 3*(-1 + a2)*a2*powf(B,2) -
				     2*powf(a,4)*B*C - 2*(-1 + a2)*a2*B*C + powf(a,4)*powf(C,2) -
				     (-1 + a2)*a2*powf(C,2) - 8*(-1 + a2)*a2*B*D -
				     4*(-1 + a2)*a2*powf(D,2) +
				     ((-1 + 2*a2)*B - a2*C)*
				     sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
					   4*(-1 + a2)*a2*powf(B + D,2))) +
		 (1 - a2)*a2*((1 - a2)*(powf(B,2)*(B + C)*(3*B + C) +
					4*powf(B,2)*(3*B + C)*D + 2*B*(7*B + C)*powf(D,2) +
					8*B*powf(D,3) + 2*powf(D,4)) +
			      (2*B + C)*(B*(B + C) + 4*B*D + 2*powf(D,2))*
			      sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
				    4*(-1 + a2)*a2*powf(B + D,2))) +
		 A*(1 - a2)*((-2 + 6*a2)*powf(B,3) + 4*a2*powf(B,2)*(-C + D) -
			     2*a2*C*(a2*powf(C,2) + powf(D,2)) +
			     2*a2*B*((-1 + 2*a2)*powf(C,2) - 2*C*D + powf(D,2)) +
			     sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
				   4*(-1 + a2)*a2*powf(B + D,2))*
			     ((-1 + 4*a2)*powf(B,2) + 4*a2*B*D -
			      a2*(powf(C,2) - 2*powf(D,2)))) +
		 powf(a,4)*(powf(B - C,2)*(B + C)*
			    sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
				  4*(-1 + a2)*a2*powf(B + D,2)) +
			    (1 - a2)*(3*powf(B,4) + 2*powf(B,3)*(C + 6*D) +
				      2*B*(C + 4*D)*(powf(C,2) + powf(D,2)) +
				      2*powf(D,2)*(2*powf(C,2) + powf(D,2)) +
				      powf(B,2)*(3*powf(C,2) + 4*C*D + 14*powf(D,2))))));

	vel = ((1 - a2)*powf(-(powf(A,2)*(-1 + a2)) + 2*A*(-1 + a2)*B + A*a2*B -
			     (-1 + a2)*powf(B,2) + a2*powf(B,2) - A*a2*C + a2*B*C + 4*a2*B*D +
			     2*a2*powf(D,2) + (A + B)*sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
							    4*(-1 + a2)*a2*powf(B + D,2)),2) +
	       a2*powf(-(A*(-1 + a2)*(B - C)) + a2*powf(B - C,2) +
		       (1 - a2)*(B*(B + C) + 4*B*D + 2*powf(D,2)) +
		       (B + C)*sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
				     4*(-1 + a2)*a2*powf(B + D,2)),2))/
	    (2.*(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) - 4*(-1 + a2)*a2*powf(B + D,2))*
	     (A - A*a2 + B + a2*C + sqrtf(powf(A*(-1 + a2) + B + a2*(-2*B + C),2) -
					  4*(-1 + a2)*a2*powf(B + D,2))));
	
	plot[ia] = sf_cmplx(sqrtf(vel)*sqrtf(fabsf(1.0f-theta))*SF_SIG(b),
			    sqrtf(vel)*sqrtf(fabsf(theta))*SF_SIG(a));
    }

    sf_complexwrite(plot,na,out);

    exit(0);
}
