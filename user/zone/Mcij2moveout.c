/* Converting Cij to moveout coefficients in 3D layered orthorhombic with possible phimuthal rotation (Sripanich and Fomel, 2016) 
These are interval parameters not effective.
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

#include <math.h>
#include <rsf.h>

int main(int argc, char* argv[])
{
	int n[3],i,j,nm;
	float o[3],d[3];
    float c33,c55,c11,c66,c12,c13,c23,c22,c44, phi, t0, mul;
    float *cc33,*cc55,*cc11,*cc66,*cc12,*cc13,*cc23,*cc22,*cc44, *pphi, *tt0;
    float *ppsi02, *ppsi20, *ppsi11, *ppsi40, *ppsi31, *ppsi22, *ppsi13, *ppsi04;
    float *tepsi02, *tepsi20, *tepsi11, *tepsi40, *tepsi31, *tepsi22, *tepsi13, *tepsi04, *tet0;
    float psi02, psi20, psi11, psi40, psi31, psi22, psi13, psi04;
    float a11, a12, a22, a1111, a1112, a1122, a1222, a2222;
    bool quarticscale, eff;
    sf_file a11o, a12o, a22o, a1111o, a1112o, a1122o, a1222o, a2222o,C33,C55,C11,C66,C12,C13,C23,C22,C44,Phi;

    sf_init (argc,argv);
    a11o = sf_output("out");
    a12o = sf_output("a12o");
    a22o = sf_output("a22o");
    a1111o = sf_output("a1111o");
    a1112o = sf_output("a1112o");
    a1122o = sf_output("a1122o");
    a1222o = sf_output("a1222o");
    a2222o = sf_output("a2222o");

    C11 = sf_input("in");
    C55 = sf_input("c55");
    C33 = sf_input("c33");
    C66 = sf_input("c66");
    C12 = sf_input("c12");
    C13 = sf_input("c13");
    C23 = sf_input("c23");
    C22 = sf_input("c22");
    C44 = sf_input("c44");
    Phi = sf_input("phi");
    
    /* get 3-D grid parameters */
    if (!sf_histint(C11,"n1",n))     sf_error("No n1= in input");
    if (!sf_histint(C11,"n2",n+1))   n[1]=1;
    if (!sf_histint(C11,"n3",n+2))   n[2]=1;
    if (!sf_histfloat(C11,"d1",d))   sf_error("No d1= in input");
    if (!sf_histfloat(C11,"d2",d+1)) d[1]=1.;
    if (!sf_histfloat(C11,"d3",d+2)) d[2]=1.;
    if (!sf_histfloat(C11,"o1",o))   o[0]=0.;
    if (!sf_histfloat(C11,"o2",o+1)) o[1]=0.;
    if (!sf_histfloat(C11,"o3",o+2)) o[2]=0.;
    
     if (!sf_getfloat("scalecij",&mul)) mul=1;
     /* Scaling of input Cij in case of GPa or km^2/s^2*/
     
     if (!sf_getbool("scalequartic",&quarticscale)) quarticscale=false;
     /* Scaling the output quartic coefficients y--multiplied by 2 t0^2 (t0 = two-way) to look at the property of the layer -> input for GMA*/

	if (!sf_getbool("eff",&eff)) eff=false;
     /* Output effective parameters instead of interval*/
     
    /* get all cij (default = orthorhombic) ---------------------------------------------------------*/
    nm = n[1]*n[2];

	/* Allocation */
	cc11 = sf_floatalloc(n[0]);
	cc22 = sf_floatalloc(n[0]);
	cc33 = sf_floatalloc(n[0]);
	cc44 = sf_floatalloc(n[0]);
	cc55 = sf_floatalloc(n[0]);
	cc66 = sf_floatalloc(n[0]);
	cc12 = sf_floatalloc(n[0]);
	cc13 = sf_floatalloc(n[0]);
	cc23 = sf_floatalloc(n[0]);
	pphi = sf_floatalloc(n[0]);
	
	tt0 = sf_floatalloc(n[0]);
	ppsi20 = sf_floatalloc(n[0]);
	ppsi11 = sf_floatalloc(n[0]);
	ppsi02 = sf_floatalloc(n[0]);
	ppsi40 = sf_floatalloc(n[0]);
	ppsi31 = sf_floatalloc(n[0]);
	ppsi22 = sf_floatalloc(n[0]);
	ppsi13 = sf_floatalloc(n[0]);
	ppsi04 = sf_floatalloc(n[0]);
	
	tet0 = sf_floatalloc(n[0]);
	tepsi20 = sf_floatalloc(n[0]);
	tepsi11 = sf_floatalloc(n[0]);
	tepsi02 = sf_floatalloc(n[0]);
	tepsi40 = sf_floatalloc(n[0]);
	tepsi31 = sf_floatalloc(n[0]);
	tepsi22 = sf_floatalloc(n[0]);
	tepsi13 = sf_floatalloc(n[0]);
	tepsi04 = sf_floatalloc(n[0]);


for (i=0;i<nm;i++) {
	
	// Read parameters
	sf_floatread(cc33,n[0],C33); sf_floatread(cc55,n[0],C55); sf_floatread(cc11,n[0],C11); 
    sf_floatread(cc12,n[0],C12); sf_floatread(cc66,n[0],C66); sf_floatread(cc13,n[0],C13); 
    sf_floatread(cc22,n[0],C22); sf_floatread(cc23,n[0],C23); sf_floatread(cc44,n[0],C44); 
    sf_floatread(pphi,n[0],Phi); 
	
	
	for (j=0;j<n[0];j++) {
	
    c11 = cc11[j]; c22 = cc22[j]; c33 = cc33[j]; c44 = cc44[j]; c55 = cc55[j];
    c66 = cc66[j]; c12 = cc12[j]; c13 = cc13[j]; c23 = cc23[j]; phi = pphi[j];
    
    c11 *=mul; c22 *=mul; c33 *=mul; c44 *=mul; c55 *=mul;
    c66 *=mul; c13 *=mul; c12 *=mul; c23 *=mul; 
    
	
	ppsi20[j] = -(pow(c33,-0.5)*(2*c23*c33*c44 + 2*c13*(c33 - c44)*c55 - 2*c23*c44*c55 - 2*c33*c44*c55 + (c33 - c44)*pow(c13,2) + c33*pow(c23,2) - c55*pow(c23,2) + c44*pow(c33,2) + c55*pow(c33,2) + cos(2*phi)*(2*c13*(c33 - c44)*c55 + 2*c23*c44*(-c33 + c55) + (c33 - c44)*pow(c13,2) + (-c33 + c55)*pow(c23,2) + (-c44 + c55)*pow(c33,2)))*pow(c33 - c44,-1)*pow(c33 - c55,-1))/2.;
	
	ppsi11[j] = cos(phi)*pow(c33,-0.5)*(2*c23*c44*(c33 - c55) + 2*c13*(-c33 + c44)*c55 + (-c33 + c44)*pow(c13,2) + (c33 - c55)*pow(c23,2) + (c44 - c55)*pow(c33,2))*pow(c33 - c44,-1)*pow(c33 - c55,-1)*sin(phi);
	
	ppsi02[j] = (pow(c33,-0.5)*(-2*c23*c33*c44 + 2*c23*c44*c55 + 2*c33*c44*c55 + 2*c13*(-c33 + c44)*c55 + (-c33 + c44)*pow(c13,2) - c33*pow(c23,2) + c55*pow(c23,2) - c44*pow(c33,2) - c55*pow(c33,2) + cos(2*phi)*(2*c13*(c33 - c44)*c55 + 2*c23*c44*(-c33 + c55) + (c33 - c44)*pow(c13,2) + (-c33 + c55)*pow(c23,2) + (-c44 + c55)*pow(c33,2)))*pow(c33 - c44,-1)*pow(c33 - c55,-1))/2.;
	
	
	ppsi40[j] = (3*pow(c33,-0.5)*pow(c33 - c44,-3)*pow(c33 - c55,-3)*(-9*c44*c55*pow(c13,4)*pow(c33,2) + 8*c44*c55*pow(c13,2)*pow(c23,2)*pow(c33,2) - 9*c44*c55*pow(c23,4)*pow(c33,2) - 
32*c12*c13*c23*c44*c55*pow(c33,3) - 32*c13*c23*c44*c55*c66*pow(c33,3) - 36*c11*c44*c55*pow(c13,2)*pow(c33,3) - 16*c23*c44*c55*pow(c13,2)*pow(c33,3) - 
12*c44*c55*c66*pow(c13,2)*pow(c33,3) - 108*c44*c55*pow(c13,3)*pow(c33,3) - 27*c44*pow(c13,4)*pow(c33,3) + 3*c55*pow(c13,4)*pow(c33,3) - 
16*c13*c44*c55*pow(c23,2)*pow(c33,3) - 36*c22*c44*c55*pow(c23,2)*pow(c33,3) - 12*c44*c55*c66*pow(c23,2)*pow(c33,3) - 8*c44*pow(c13,2)*pow(c23,2)*pow(c33,3) - 8*c55*pow(c13,2)*pow(c23,2)*pow(c33,3) - 108*c44*c55*pow(c23,3)*pow(c33,3) + 3*c44*pow(c23,4)*pow(c33,3) - 27*c55*pow(c23,4)*pow(c33,3) + 
16*c12*c13*c23*c44*pow(c33,4) + 16*c12*c13*c23*c55*pow(c33,4) + 72*c11*c13*c44*c55*pow(c33,4) + 16*c12*c13*c44*c55*pow(c33,4) + 16*c12*c23*c44*c55*pow(c33,4) + 
24*c13*c23*c44*c55*pow(c33,4) + 72*c22*c23*c44*c55*pow(c33,4) + 16*c13*c23*c44*c66*pow(c33,4) + 16*c13*c23*c55*c66*pow(c33,4) + 40*c13*c44*c55*c66*pow(c33,4) + 
40*c23*c44*c55*c66*pow(c33,4) + 36*c11*c44*pow(c13,2)*pow(c33,4) + 12*c23*c44*pow(c13,2)*pow(c33,4) + 12*c11*c55*pow(c13,2)*pow(c33,4) - 
18*c44*c55*pow(c13,2)*pow(c33,4) + 12*c44*c66*pow(c13,2)*pow(c33,4) + 4*c55*c66*pow(c13,2)*pow(c33,4) + 36*c55*pow(c13,3)*pow(c33,4) + 9*pow(c13,4)*pow(c33,4) + 
12*c22*c44*pow(c23,2)*pow(c33,4) + 12*c13*c55*pow(c23,2)*pow(c33,4) + 36*c22*c55*pow(c23,2)*pow(c33,4) - 18*c44*c55*pow(c23,2)*pow(c33,4) + 
4*c44*c66*pow(c23,2)*pow(c33,4) + 12*c55*c66*pow(c23,2)*pow(c33,4) + 6*pow(c13,2)*pow(c23,2)*pow(c33,4) + 36*c44*pow(c23,3)*pow(c33,4) + 9*pow(c23,4)*pow(c33,4) - 
8*c12*c13*c23*pow(c33,5) - 8*c12*c13*c44*pow(c33,5) - 24*c22*c23*c44*pow(c33,5) - 24*c11*c13*c55*pow(c33,5) - 8*c12*c23*c55*pow(c33,5) - 
8*c12*c44*c55*pow(c33,5) + 4*c13*c44*c55*pow(c33,5) + 4*c23*c44*c55*pow(c33,5) - 8*c13*c23*c66*pow(c33,5) - 8*c13*c44*c66*pow(c33,5) - 8*c23*c44*c66*pow(c33,5) - 
8*c13*c55*c66*pow(c33,5) - 8*c23*c55*c66*pow(c33,5) - 8*c44*c55*c66*pow(c33,5) - 12*c11*pow(c13,2)*pow(c33,5) + 2*c44*pow(c13,2)*pow(c33,5) + 
6*c55*pow(c13,2)*pow(c33,5) - 4*c66*pow(c13,2)*pow(c33,5) - 12*c22*pow(c23,2)*pow(c33,5) + 6*c44*pow(c23,2)*pow(c33,5) + 2*c55*pow(c23,2)*pow(c33,5) - 
4*c66*pow(c23,2)*pow(c33,5) - 2*c44*c55*pow(c33,6) + 9*c33*c55*pow(c13,4)*pow(c44,2) + 16*c12*c13*c23*c55*pow(c33,2)*pow(c44,2) + 
16*c13*c23*c55*c66*pow(c33,2)*pow(c44,2) + 36*c11*c55*pow(c13,2)*pow(c33,2)*pow(c44,2) + 16*c23*c55*pow(c13,2)*pow(c33,2)*pow(c44,2) + 
12*c55*c66*pow(c13,2)*pow(c33,2)*pow(c44,2) + 108*c55*pow(c13,3)*pow(c33,2)*pow(c44,2) + 27*pow(c13,4)*pow(c33,2)*pow(c44,2) + 
4*c13*c55*pow(c23,2)*pow(c33,2)*pow(c44,2) + 2*pow(c13,2)*pow(c23,2)*pow(c33,2)*pow(c44,2) - 36*c55*pow(c23,3)*pow(c33,2)*pow(c44,2) - 
8*c12*c13*c23*pow(c33,3)*pow(c44,2) - 72*c11*c13*c55*pow(c33,3)*pow(c44,2) - 32*c12*c13*c55*pow(c33,3)*pow(c44,2) - 8*c12*c23*c55*pow(c33,3)*pow(c44,2) - 
32*c13*c23*c55*pow(c33,3)*pow(c44,2) - 72*c22*c23*c55*pow(c33,3)*pow(c44,2) - 8*c13*c23*c66*pow(c33,3)*pow(c44,2) - 56*c13*c55*c66*pow(c33,3)*pow(c44,2) - 
32*c23*c55*c66*pow(c33,3)*pow(c44,2) - 36*c11*pow(c13,2)*pow(c33,3)*pow(c44,2) - 16*c23*pow(c13,2)*pow(c33,3)*pow(c44,2) + 
10*c55*pow(c13,2)*pow(c33,3)*pow(c44,2) - 12*c66*pow(c13,2)*pow(c33,3)*pow(c44,2) - 164*c55*pow(c23,2)*pow(c33,3)*pow(c44,2) + 
12*pow(c23,3)*pow(c33,3)*pow(c44,2) + 16*c12*c13*pow(c33,4)*pow(c44,2) + 24*c22*c23*pow(c33,4)*pow(c44,2) + 16*c12*c55*pow(c33,4)*pow(c44,2) + 
36*c22*c55*pow(c33,4)*pow(c44,2) - 36*c23*c55*pow(c33,4)*pow(c44,2) + 16*c13*c66*pow(c33,4)*pow(c44,2) + 8*c23*c66*pow(c33,4)*pow(c44,2) + 
28*c55*c66*pow(c33,4)*pow(c44,2) + 54*pow(c23,2)*pow(c33,4)*pow(c44,2) - 12*c22*pow(c33,5)*pow(c44,2) + 12*c23*pow(c33,5)*pow(c44,2) + 
17*c55*pow(c33,5)*pow(c44,2) - 4*c66*pow(c33,5)*pow(c44,2) - 3*pow(c33,6)*pow(c44,2) - 12*c11*c33*c55*pow(c13,2)*pow(c44,3) - 
4*c33*c55*c66*pow(c13,2)*pow(c44,3) - 36*c33*c55*pow(c13,3)*pow(c44,3) - 9*c33*pow(c13,4)*pow(c44,3) - 3*c55*pow(c13,4)*pow(c44,3) + 
24*c11*c13*c55*pow(c33,2)*pow(c44,3) + 16*c12*c13*c55*pow(c33,2)*pow(c44,3) + 8*c13*c23*c55*pow(c33,2)*pow(c44,3) + 24*c13*c55*c66*pow(c33,2)*pow(c44,3) + 
12*c11*pow(c13,2)*pow(c33,2)*pow(c44,3) + 4*c23*pow(c13,2)*pow(c33,2)*pow(c44,3) + 2*c55*pow(c13,2)*pow(c33,2)*pow(c44,3) + 
4*c66*pow(c13,2)*pow(c33,2)*pow(c44,3) - 36*c55*pow(c23,2)*pow(c33,2)*pow(c44,3) - 8*c12*c13*pow(c33,3)*pow(c44,3) - 8*c12*c55*pow(c33,3)*pow(c44,3) - 
4*c13*c55*pow(c33,3)*pow(c44,3) - 36*c22*c55*pow(c33,3)*pow(c44,3) - 112*c23*c55*pow(c33,3)*pow(c44,3) - 8*c13*c66*pow(c33,3)*pow(c44,3) - 
20*c55*c66*pow(c33,3)*pow(c44,3) - 2*pow(c13,2)*pow(c33,3)*pow(c44,3) + 12*pow(c23,2)*pow(c33,3)*pow(c44,3) + 12*c22*pow(c33,4)*pow(c44,3) + 
36*c23*pow(c33,4)*pow(c44,3) - 51*c55*pow(c33,4)*pow(c44,3) + 4*c66*pow(c33,4)*pow(c44,3) + 15*pow(c33,5)*pow(c44,3) + 9*c33*c44*pow(c23,4)*pow(c55,2) + 
16*c12*c13*c23*c44*pow(c33,2)*pow(c55,2) + 16*c13*c23*c44*c66*pow(c33,2)*pow(c55,2) + 4*c23*c44*pow(c13,2)*pow(c33,2)*pow(c55,2) - 
36*c44*pow(c13,3)*pow(c33,2)*pow(c55,2) + 16*c13*c44*pow(c23,2)*pow(c33,2)*pow(c55,2) + 36*c22*c44*pow(c23,2)*pow(c33,2)*pow(c55,2) + 
12*c44*c66*pow(c23,2)*pow(c33,2)*pow(c55,2) + 2*pow(c13,2)*pow(c23,2)*pow(c33,2)*pow(c55,2) + 108*c44*pow(c23,3)*pow(c33,2)*pow(c55,2) + 
27*pow(c23,4)*pow(c33,2)*pow(c55,2) - 8*c12*c13*c23*pow(c33,3)*pow(c55,2) - 72*c11*c13*c44*pow(c33,3)*pow(c55,2) - 8*c12*c13*c44*pow(c33,3)*pow(c55,2) - 
32*c12*c23*c44*pow(c33,3)*pow(c55,2) - 32*c13*c23*c44*pow(c33,3)*pow(c55,2) - 72*c22*c23*c44*pow(c33,3)*pow(c55,2) - 8*c13*c23*c66*pow(c33,3)*pow(c55,2) - 
32*c13*c44*c66*pow(c33,3)*pow(c55,2) - 56*c23*c44*c66*pow(c33,3)*pow(c55,2) - 164*c44*pow(c13,2)*pow(c33,3)*pow(c55,2) + 12*pow(c13,3)*pow(c33,3)*pow(c55,2) - 16*c13*pow(c23,2)*pow(c33,3)*pow(c55,2) - 36*c22*pow(c23,2)*pow(c33,3)*pow(c55,2) + 10*c44*pow(c23,2)*pow(c33,3)*pow(c55,2) - 
12*c66*pow(c23,2)*pow(c33,3)*pow(c55,2) + 24*c11*c13*pow(c33,4)*pow(c55,2) + 16*c12*c23*pow(c33,4)*pow(c55,2) + 36*c11*c44*pow(c33,4)*pow(c55,2) + 
16*c12*c44*pow(c33,4)*pow(c55,2) - 36*c13*c44*pow(c33,4)*pow(c55,2) + 8*c13*c66*pow(c33,4)*pow(c55,2) + 16*c23*c66*pow(c33,4)*pow(c55,2) + 
28*c44*c66*pow(c33,4)*pow(c55,2) + 54*pow(c13,2)*pow(c33,4)*pow(c55,2) - 12*c11*pow(c33,5)*pow(c55,2) + 12*c13*pow(c33,5)*pow(c55,2) + 
17*c44*pow(c33,5)*pow(c55,2) - 4*c66*pow(c33,5)*pow(c55,2) - 3*pow(c33,6)*pow(c55,2) - 8*c12*c13*c23*c33*pow(c44,2)*pow(c55,2) - 
8*c13*c23*c33*c66*pow(c44,2)*pow(c55,2) + 36*c33*pow(c13,3)*pow(c44,2)*pow(c55,2) - 2*pow(c13,2)*pow(c23,2)*pow(c44,2)*pow(c55,2) + 
36*c33*pow(c23,3)*pow(c44,2)*pow(c55,2) + 72*c11*c13*pow(c33,2)*pow(c44,2)*pow(c55,2) + 16*c12*c13*pow(c33,2)*pow(c44,2)*pow(c55,2) + 
16*c12*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 32*c13*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 72*c22*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 
40*c13*c66*pow(c33,2)*pow(c44,2)*pow(c55,2) + 40*c23*c66*pow(c33,2)*pow(c44,2)*pow(c55,2) + 170*pow(c13,2)*pow(c33,2)*pow(c44,2)*pow(c55,2) + 
170*pow(c23,2)*pow(c33,2)*pow(c44,2)*pow(c55,2) - 36*c11*pow(c33,3)*pow(c44,2)*pow(c55,2) - 32*c12*pow(c33,3)*pow(c44,2)*pow(c55,2) + 
20*c13*pow(c33,3)*pow(c44,2)*pow(c55,2) - 36*c22*pow(c33,3)*pow(c44,2)*pow(c55,2) + 20*c23*pow(c33,3)*pow(c44,2)*pow(c55,2) - 
56*c66*pow(c33,3)*pow(c44,2)*pow(c55,2) - 42*pow(c33,4)*pow(c44,2)*pow(c55,2) - 24*c11*c13*c33*pow(c44,3)*pow(c55,2) - 8*c12*c13*c33*pow(c44,3)*pow(c55,2) - 16*c13*c33*c66*pow(c44,3)*pow(c55,2) - 4*c23*pow(c13,2)*pow(c44,3)*pow(c55,2) - 60*c33*pow(c13,2)*pow(c44,3)*pow(c55,2) - 12*pow(c13,3)*pow(c44,3)*pow(c55,2) + 
36*c33*pow(c23,2)*pow(c44,3)*pow(c55,2) + 12*c11*pow(c33,2)*pow(c44,3)*pow(c55,2) + 16*c12*pow(c33,2)*pow(c44,3)*pow(c55,2) + 
4*c13*pow(c33,2)*pow(c44,3)*pow(c55,2) + 36*c22*pow(c33,2)*pow(c44,3)*pow(c55,2) + 124*c23*pow(c33,2)*pow(c44,3)*pow(c55,2) + 
32*c66*pow(c33,2)*pow(c44,3)*pow(c55,2) + 64*pow(c33,3)*pow(c44,3)*pow(c55,2) - 12*c22*c33*c44*pow(c23,2)*pow(c55,3) - 4*c33*c44*c66*pow(c23,2)*pow(c55,3) - 36*c33*c44*pow(c23,3)*pow(c55,3) - 9*c33*pow(c23,4)*pow(c55,3) - 3*c44*pow(c23,4)*pow(c55,3) + 16*c12*c23*c44*pow(c33,2)*pow(c55,3) + 
8*c13*c23*c44*pow(c33,2)*pow(c55,3) + 24*c22*c23*c44*pow(c33,2)*pow(c55,3) + 24*c23*c44*c66*pow(c33,2)*pow(c55,3) - 36*c44*pow(c13,2)*pow(c33,2)*pow(c55,3) + 4*c13*pow(c23,2)*pow(c33,2)*pow(c55,3) + 12*c22*pow(c23,2)*pow(c33,2)*pow(c55,3) + 2*c44*pow(c23,2)*pow(c33,2)*pow(c55,3) + 
4*c66*pow(c23,2)*pow(c33,2)*pow(c55,3) - 8*c12*c23*pow(c33,3)*pow(c55,3) - 36*c11*c44*pow(c33,3)*pow(c55,3) - 8*c12*c44*pow(c33,3)*pow(c55,3) - 
112*c13*c44*pow(c33,3)*pow(c55,3) - 4*c23*c44*pow(c33,3)*pow(c55,3) - 8*c23*c66*pow(c33,3)*pow(c55,3) - 20*c44*c66*pow(c33,3)*pow(c55,3) + 
12*pow(c13,2)*pow(c33,3)*pow(c55,3) - 2*pow(c23,2)*pow(c33,3)*pow(c55,3) + 12*c11*pow(c33,4)*pow(c55,3) + 36*c13*pow(c33,4)*pow(c55,3) - 
51*c44*pow(c33,4)*pow(c55,3) + 4*c66*pow(c33,4)*pow(c55,3) + 15*pow(c33,5)*pow(c55,3) - 8*c12*c23*c33*pow(c44,2)*pow(c55,3) - 
24*c22*c23*c33*pow(c44,2)*pow(c55,3) - 16*c23*c33*c66*pow(c44,2)*pow(c55,3) + 36*c33*pow(c13,2)*pow(c44,2)*pow(c55,3) - 4*c13*pow(c23,2)*pow(c44,2)*pow(c55,3) - 
60*c33*pow(c23,2)*pow(c44,2)*pow(c55,3) - 12*pow(c23,3)*pow(c44,2)*pow(c55,3) + 36*c11*pow(c33,2)*pow(c44,2)*pow(c55,3) + 
16*c12*pow(c33,2)*pow(c44,2)*pow(c55,3) + 124*c13*pow(c33,2)*pow(c44,2)*pow(c55,3) + 12*c22*pow(c33,2)*pow(c44,2)*pow(c55,3) + 
4*c23*pow(c33,2)*pow(c44,2)*pow(c55,3) + 32*c66*pow(c33,2)*pow(c44,2)*pow(c55,3) + 64*pow(c33,3)*pow(c44,2)*pow(c55,3) - 8*c13*c23*pow(c44,3)*pow(c55,3) - 
12*c11*c33*pow(c44,3)*pow(c55,3) - 8*c12*c33*pow(c44,3)*pow(c55,3) - 48*c13*c33*pow(c44,3)*pow(c55,3) - 12*c22*c33*pow(c44,3)*pow(c55,3) - 
48*c23*c33*pow(c44,3)*pow(c55,3) - 16*c33*c66*pow(c44,3)*pow(c55,3) - 12*pow(c13,2)*pow(c44,3)*pow(c55,3) - 12*pow(c23,2)*pow(c44,3)*pow(c55,3) - 
40*pow(c33,2)*pow(c44,3)*pow(c55,3) + 4*cos(2*phi)*(4*c13*c33*c55*(-2*c11*(c33 - c55) + c55*(c33 + 3*c55))*pow(c33 - c44,3) + 
4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + (3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) - 
4*c23*c33*c44*(-2*c22*(c33 - c44) + c44*(c33 + 3*c44))*pow(c33 - c55,3) - 4*c44*(3*c33 + c44)*pow(c23,3)*pow(c33 - c55,3) - 
(3*c33 + c44)*pow(c23,4)*pow(c33 - c55,3) + 2*pow(c23,2)*(2*c22*c33*(c33 - c44) - c44*(9*c33*c44 + pow(c33,2) + 2*pow(c44,2)))*pow(c33 - c55,3) - 
2*pow(c13,2)*pow(c33 - c44,3)*(2*c11*c33*(c33 - c55) - c55*(9*c33*c55 + pow(c33,2) + 2*pow(c55,2))) + 
c33*(4*c22*(c33 - c44)*pow(c44,2)*pow(c33 - c55,3) + pow(c33,5)*(pow(c44,2) - pow(c55,2)) - 
2*c44*(7*c44*(c44 - c55) + 6*c11*(c44 + c55))*pow(c33,2)*pow(c55,2) + 4*c11*c33*(c44 + 3*c55)*pow(c44,2)*pow(c55,2) + 
c55*pow(c33,3)*(3*c44*(4*c11 - 5*c55)*c55 + 15*pow(c44,3) + 4*c11*pow(c55,2)) + 
pow(c33,4)*(-3*c55*pow(c44,2) - 5*pow(c44,3) + 3*c44*pow(c55,2) + (-4*c11 + 5*c55)*pow(c55,2)) - 4*c11*pow(c44,3)*pow(c55,3))) + 
cos(4*phi)*(4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + (3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) + 4*c44*(3*c33 + c44)*pow(c23,3)*pow(c33 - c55,3) + 
(3*c33 + c44)*pow(c23,4)*pow(c33 - c55,3) - 2*(c33 - c55)*pow(c23,2)*
((-c44 + c55 - 2*c66)*pow(c33,4) + pow(c33,3)*(2*c44*(c55 + c66) + c55*(c55 + 4*c66) - 9*pow(c44,2)) + 2*c33*c44*c55*(-3*c44*c55 + c55*c66 + 2*pow(c44,2)) + 2*c22*c33*(c33 - c44)*pow(c33 - c55,2) - 2*pow(c44,3)*pow(c55,2) - pow(c33,2)*(c44*c55*(5*c55 + 4*c66) - 17*c55*pow(c44,2) + 2*pow(c44,3) + 2*c66*pow(c55,2))
) - 4*c23*c33*(c33 - c55)*(-(c33*c44*c55*(5*c44*c55 + 2*c12*(c44 + 2*c55) + 6*c44*c66 + 6*c55*c66 - 5*pow(c44,2))) - 
pow(c33,3)*(-(c44*(c55 - 2*c66)) + 2*c55*(c12 + c66) + pow(c44,2)) + 2*c22*(c33 - c44)*c44*pow(c33 - c55,2) + 2*(c12 + 2*c66)*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(c44*c55*(4*c12 + c55 + 8*c66) + 2*(c55 + c66)*pow(c44,2) - 3*pow(c44,3) + 2*(c12 + c66)*pow(c55,2))) + 
4*c13*(c33 - c44)*(-((c33 - c55)*c55*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2))) + 2*c12*c33*(c33 - c44)*(c23 + c44)*pow(c33 - c55,2) + 
2*c23*(c33 - c55)*(c33*c44*c55*(c44 + c55 + c66) - (c55*c66 + c44*(3*c55 + c66))*pow(c33,2) + c66*pow(c33,3) + pow(c44,2)*pow(c55,2)) - 
c33*((c44*(c55 - 2*c66) - c55*(c55 + 2*c66))*pow(c33,3) + 2*c11*(c33 - c55)*c55*pow(c33 - c44,2) - 
c33*c44*c55*(5*c44*c55 + 6*c44*c66 + 6*c55*c66 - 5*pow(c55,2)) + 4*c66*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(2*c44*c55*(c55 + 4*c66) + (c55 + 2*c66)*pow(c44,2) + (-3*c55 + 2*c66)*pow(c55,2)))) - 
2*(c33 - c44)*pow(c13,2)*(-2*c23*c44*(c33 - c55)*(c44*c55 + c33*(c44 + c55) - 3*pow(c33,2)) - 4*c44*c55*c66*pow(c33,2) + 
(c33 - c55)*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2)) + 2*c44*c55*pow(c33,3) + 4*c44*c66*pow(c33,3) + 2*c55*c66*pow(c33,3) + c44*pow(c33,4) - c55*pow(c33,4) - 2*c66*pow(c33,4) + 2*c11*c33*(c33 - c55)*pow(c33 - c44,2) + 2*c33*c55*c66*pow(c44,2) - 5*c55*pow(c33,2)*pow(c44,2) - 
2*c66*pow(c33,2)*pow(c44,2) + pow(c33,3)*pow(c44,2) + 17*c44*pow(c33,2)*pow(c55,2) - 9*pow(c33,3)*pow(c55,2) - 6*c33*pow(c44,2)*pow(c55,2) + 
4*c33*c44*pow(c55,3) - 2*pow(c33,2)*pow(c55,3) - 2*pow(c44,2)*pow(c55,3)) - 
c33*(4*c22*(c33 - c44)*pow(c44,2)*pow(c33 - c55,3) + pow(c33,5)*pow(c44 - c55,2) - 
4*c33*(c11*(c44 + 3*c55) - 4*(c44 + c55)*(c12 + 2*c66))*pow(c44,2)*pow(c55,2) - 
pow(c33,4)*(c44*c55*(8*c12 - 5*c55 + 8*c66) + (-5*c55 + 4*c66)*pow(c44,2) + 5*pow(c44,3) + (-4*c11 + 5*c55 + 4*c66)*pow(c55,2)) - 
4*c44*c55*pow(c33,2)*(-3*c11*c55*(c44 + c55) + 2*c12*(4*c44*c55 + pow(c44,2) + pow(c55,2)) + c66*(14*c44*c55 + 5*pow(c44,2) + 5*pow(c55,2))) + 
4*(c11 - 2*(c12 + 2*c66))*pow(c44,3)*pow(c55,3) + pow(c33,3)*
(16*c12*c44*c55*(c44 + c55) + 2*c55*(-9*c55 + 14*c66)*pow(c44,2) + (9*c55 + 4*c66)*pow(c44,3) + c44*(-12*c11 + 9*c55 + 28*c66)*pow(c55,2) + 
4*(-c11 + c66)*pow(c55,3))))))/8.;
	
	
	
	ppsi31[j] = (3*pow(c33,-0.5)*pow(c33 - c44,-3)*pow(c33 - c55,-3)*(3*c44*c55*pow(c23,4)*pow(c33,2) + 12*c22*c44*c55*pow(c23,2)*pow(c33,3) + 36*c44*c55*pow(c23,3)*pow(c33,3) - 
c44*pow(c23,4)*pow(c33,3) + 9*c55*pow(c23,4)*pow(c33,3) - 24*c22*c23*c44*c55*pow(c33,4) - 4*c22*c44*pow(c23,2)*pow(c33,4) - 12*c22*c55*pow(c23,2)*pow(c33,4) + 6*c44*c55*pow(c23,2)*pow(c33,4) - 12*c44*pow(c23,3)*pow(c33,4) - 3*pow(c23,4)*pow(c33,4) + 8*c22*c23*c44*pow(c33,5) + 4*c22*pow(c23,2)*pow(c33,5) - 
2*c44*pow(c23,2)*pow(c33,5) + 4*c13*c33*c55*(-2*c11*(c33 - c55) + c55*(c33 + 3*c55))*pow(c33 - c44,3) + 4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + 
(3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) + 12*c55*pow(c23,3)*pow(c33,2)*pow(c44,2) + 24*c22*c23*c55*pow(c33,3)*pow(c44,2) + 
54*c55*pow(c23,2)*pow(c33,3)*pow(c44,2) - 4*pow(c23,3)*pow(c33,3)*pow(c44,2) - 8*c22*c23*pow(c33,4)*pow(c44,2) - 12*c22*c55*pow(c33,4)*pow(c44,2) + 
12*c23*c55*pow(c33,4)*pow(c44,2) - 18*pow(c23,2)*pow(c33,4)*pow(c44,2) + 4*c22*pow(c33,5)*pow(c44,2) - 4*c23*pow(c33,5)*pow(c44,2) - 3*c55*pow(c33,5)*pow(c44,2) + 
pow(c33,6)*pow(c44,2) + 12*c55*pow(c23,2)*pow(c33,2)*pow(c44,3) + 12*c22*c55*pow(c33,3)*pow(c44,3) + 36*c23*c55*pow(c33,3)*pow(c44,3) - 
4*pow(c23,2)*pow(c33,3)*pow(c44,3) - 4*c22*pow(c33,4)*pow(c44,3) - 12*c23*pow(c33,4)*pow(c44,3) + 15*c55*pow(c33,4)*pow(c44,3) - 5*pow(c33,5)*pow(c44,3) - 
3*c33*c44*pow(c23,4)*pow(c55,2) - 12*c22*c44*pow(c23,2)*pow(c33,2)*pow(c55,2) - 36*c44*pow(c23,3)*pow(c33,2)*pow(c55,2) - 9*pow(c23,4)*pow(c33,2)*pow(c55,2) + 24*c22*c23*c44*pow(c33,3)*pow(c55,2) + 12*c22*pow(c23,2)*pow(c33,3)*pow(c55,2) - 6*c44*pow(c23,2)*pow(c33,3)*pow(c55,2) + 12*c11*c44*pow(c33,4)*pow(c55,2) - 4*c11*pow(c33,5)*pow(c55,2) + 3*c44*pow(c33,5)*pow(c55,2) - pow(c33,6)*pow(c55,2) - 12*c33*pow(c23,3)*pow(c44,2)*pow(c55,2) - 
24*c22*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) - 54*pow(c23,2)*pow(c33,2)*pow(c44,2)*pow(c55,2) - 12*c11*pow(c33,3)*pow(c44,2)*pow(c55,2) + 
12*c22*pow(c33,3)*pow(c44,2)*pow(c55,2) - 12*c23*pow(c33,3)*pow(c44,2)*pow(c55,2) - 12*c33*pow(c23,2)*pow(c44,3)*pow(c55,2) + 
4*c11*pow(c33,2)*pow(c44,3)*pow(c55,2) - 12*c22*pow(c33,2)*pow(c44,3)*pow(c55,2) - 36*c23*pow(c33,2)*pow(c44,3)*pow(c55,2) - 14*pow(c33,3)*pow(c44,3)*pow(c55,2) - 
2*pow(c13,2)*pow(c33 - c44,3)*(2*c11*c33*(c33 - c55) - c55*(9*c33*c55 + pow(c33,2) + 2*pow(c55,2))) + 4*c22*c33*c44*pow(c23,2)*pow(c55,3) + 
12*c33*c44*pow(c23,3)*pow(c55,3) + 3*c33*pow(c23,4)*pow(c55,3) + c44*pow(c23,4)*pow(c55,3) - 8*c22*c23*c44*pow(c33,2)*pow(c55,3) - 
4*c22*pow(c23,2)*pow(c33,2)*pow(c55,3) + 2*c44*pow(c23,2)*pow(c33,2)*pow(c55,3) - 12*c11*c44*pow(c33,3)*pow(c55,3) + 4*c11*pow(c33,4)*pow(c55,3) - 
15*c44*pow(c33,4)*pow(c55,3) + 5*pow(c33,5)*pow(c55,3) + 8*c22*c23*c33*pow(c44,2)*pow(c55,3) + 18*c33*pow(c23,2)*pow(c44,2)*pow(c55,3) + 
4*pow(c23,3)*pow(c44,2)*pow(c55,3) + 12*c11*pow(c33,2)*pow(c44,2)*pow(c55,3) - 4*c22*pow(c33,2)*pow(c44,2)*pow(c55,3) + 4*c23*pow(c33,2)*pow(c44,2)*pow(c55,3) + 
14*pow(c33,3)*pow(c44,2)*pow(c55,3) - 4*c11*c33*pow(c44,3)*pow(c55,3) + 4*c22*c33*pow(c44,3)*pow(c55,3) + 12*c23*c33*pow(c44,3)*pow(c55,3) + 
4*pow(c23,2)*pow(c44,3)*pow(c55,3) + cos(2*phi)*(4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + (3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) + 
4*c44*(3*c33 + c44)*pow(c23,3)*pow(c33 - c55,3) + (3*c33 + c44)*pow(c23,4)*pow(c33 - c55,3) - 
2*(c33 - c55)*pow(c23,2)*((-c44 + c55 - 2*c66)*pow(c33,4) + pow(c33,3)*(2*c44*(c55 + c66) + c55*(c55 + 4*c66) - 9*pow(c44,2)) + 
2*c33*c44*c55*(-3*c44*c55 + c55*c66 + 2*pow(c44,2)) + 2*c22*c33*(c33 - c44)*pow(c33 - c55,2) - 2*pow(c44,3)*pow(c55,2) - 
pow(c33,2)*(c44*c55*(5*c55 + 4*c66) - 17*c55*pow(c44,2) + 2*pow(c44,3) + 2*c66*pow(c55,2))) - 
4*c23*c33*(c33 - c55)*(-(c33*c44*c55*(5*c44*c55 + 2*c12*(c44 + 2*c55) + 6*c44*c66 + 6*c55*c66 - 5*pow(c44,2))) - 
pow(c33,3)*(-(c44*(c55 - 2*c66)) + 2*c55*(c12 + c66) + pow(c44,2)) + 2*c22*(c33 - c44)*c44*pow(c33 - c55,2) + 2*(c12 + 2*c66)*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(c44*c55*(4*c12 + c55 + 8*c66) + 2*(c55 + c66)*pow(c44,2) - 3*pow(c44,3) + 2*(c12 + c66)*pow(c55,2))) + 
4*c13*(c33 - c44)*(-((c33 - c55)*c55*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2))) + 2*c12*c33*(c33 - c44)*(c23 + c44)*pow(c33 - c55,2) + 
2*c23*(c33 - c55)*(c33*c44*c55*(c44 + c55 + c66) - (c55*c66 + c44*(3*c55 + c66))*pow(c33,2) + c66*pow(c33,3) + pow(c44,2)*pow(c55,2)) - 
c33*((c44*(c55 - 2*c66) - c55*(c55 + 2*c66))*pow(c33,3) + 2*c11*(c33 - c55)*c55*pow(c33 - c44,2) - 
c33*c44*c55*(5*c44*c55 + 6*c44*c66 + 6*c55*c66 - 5*pow(c55,2)) + 4*c66*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(2*c44*c55*(c55 + 4*c66) + (c55 + 2*c66)*pow(c44,2) + (-3*c55 + 2*c66)*pow(c55,2)))) - 
2*(c33 - c44)*pow(c13,2)*(-2*c23*c44*(c33 - c55)*(c44*c55 + c33*(c44 + c55) - 3*pow(c33,2)) - 4*c44*c55*c66*pow(c33,2) + 
(c33 - c55)*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2)) + 2*c44*c55*pow(c33,3) + 4*c44*c66*pow(c33,3) + 2*c55*c66*pow(c33,3) + c44*pow(c33,4) - c55*pow(c33,4) - 2*c66*pow(c33,4) + 2*c11*c33*(c33 - c55)*pow(c33 - c44,2) + 2*c33*c55*c66*pow(c44,2) - 5*c55*pow(c33,2)*pow(c44,2) - 
2*c66*pow(c33,2)*pow(c44,2) + pow(c33,3)*pow(c44,2) + 17*c44*pow(c33,2)*pow(c55,2) - 9*pow(c33,3)*pow(c55,2) - 6*c33*pow(c44,2)*pow(c55,2) + 
4*c33*c44*pow(c55,3) - 2*pow(c33,2)*pow(c55,3) - 2*pow(c44,2)*pow(c55,3)) - 
c33*(4*c22*(c33 - c44)*pow(c44,2)*pow(c33 - c55,3) + pow(c33,5)*pow(c44 - c55,2) - 
4*c33*(c11*(c44 + 3*c55) - 4*(c44 + c55)*(c12 + 2*c66))*pow(c44,2)*pow(c55,2) - 
pow(c33,4)*(c44*c55*(8*c12 - 5*c55 + 8*c66) + (-5*c55 + 4*c66)*pow(c44,2) + 5*pow(c44,3) + (-4*c11 + 5*c55 + 4*c66)*pow(c55,2)) - 
4*c44*c55*pow(c33,2)*(-3*c11*c55*(c44 + c55) + 2*c12*(4*c44*c55 + pow(c44,2) + pow(c55,2)) + c66*(14*c44*c55 + 5*pow(c44,2) + 5*pow(c55,2))) + 
4*(c11 - 2*(c12 + 2*c66))*pow(c44,3)*pow(c55,3) + pow(c33,3)*
(16*c12*c44*c55*(c44 + c55) + 2*c55*(-9*c55 + 14*c66)*pow(c44,2) + (9*c55 + 4*c66)*pow(c44,3) + c44*(-12*c11 + 9*c55 + 28*c66)*pow(c55,2) + 
4*(-c11 + c66)*pow(c55,3)))))*sin(2*phi))/4.;



	ppsi22[j] = (pow(c33,-0.5)*pow(c33 - c44,-3)*pow(c33 - c55,-3)*(-9*c44*c55*pow(c23,4)*pow(c33,2) - 36*c22*c44*c55*pow(c23,2)*pow(c33,3) - 12*c44*c55*c66*pow(c23,2)*pow(c33,3) - 
108*c44*c55*pow(c23,3)*pow(c33,3) + 3*c44*pow(c23,4)*pow(c33,3) - 27*c55*pow(c23,4)*pow(c33,3) + 16*c12*c23*c44*c55*pow(c33,4) + 72*c22*c23*c44*c55*pow(c33,4) + 
40*c23*c44*c55*c66*pow(c33,4) + 12*c22*c44*pow(c23,2)*pow(c33,4) + 36*c22*c55*pow(c23,2)*pow(c33,4) - 18*c44*c55*pow(c23,2)*pow(c33,4) + 
4*c44*c66*pow(c23,2)*pow(c33,4) + 12*c55*c66*pow(c23,2)*pow(c33,4) + 36*c44*pow(c23,3)*pow(c33,4) + 9*pow(c23,4)*pow(c33,4) - 24*c22*c23*c44*pow(c33,5) - 
8*c12*c23*c55*pow(c33,5) - 8*c12*c44*c55*pow(c33,5) + 4*c23*c44*c55*pow(c33,5) - 8*c23*c44*c66*pow(c33,5) - 8*c23*c55*c66*pow(c33,5) - 8*c44*c55*c66*pow(c33,5) - 
12*c22*pow(c23,2)*pow(c33,5) + 6*c44*pow(c23,2)*pow(c33,5) + 2*c55*pow(c23,2)*pow(c33,5) - 4*c66*pow(c23,2)*pow(c33,5) - 2*c44*c55*pow(c33,6) + 
12*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + 3*(3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) - 36*c55*pow(c23,3)*pow(c33,2)*pow(c44,2) - 
8*c12*c23*c55*pow(c33,3)*pow(c44,2) - 72*c22*c23*c55*pow(c33,3)*pow(c44,2) - 32*c23*c55*c66*pow(c33,3)*pow(c44,2) - 164*c55*pow(c23,2)*pow(c33,3)*pow(c44,2) + 12*pow(c23,3)*pow(c33,3)*pow(c44,2) + 24*c22*c23*pow(c33,4)*pow(c44,2) + 16*c12*c55*pow(c33,4)*pow(c44,2) + 36*c22*c55*pow(c33,4)*pow(c44,2) - 
36*c23*c55*pow(c33,4)*pow(c44,2) + 8*c23*c66*pow(c33,4)*pow(c44,2) + 28*c55*c66*pow(c33,4)*pow(c44,2) + 54*pow(c23,2)*pow(c33,4)*pow(c44,2) - 
12*c22*pow(c33,5)*pow(c44,2) + 12*c23*pow(c33,5)*pow(c44,2) + 17*c55*pow(c33,5)*pow(c44,2) - 4*c66*pow(c33,5)*pow(c44,2) - 3*pow(c33,6)*pow(c44,2) - 
36*c55*pow(c23,2)*pow(c33,2)*pow(c44,3) - 8*c12*c55*pow(c33,3)*pow(c44,3) - 36*c22*c55*pow(c33,3)*pow(c44,3) - 112*c23*c55*pow(c33,3)*pow(c44,3) - 
20*c55*c66*pow(c33,3)*pow(c44,3) + 12*pow(c23,2)*pow(c33,3)*pow(c44,3) + 12*c22*pow(c33,4)*pow(c44,3) + 36*c23*pow(c33,4)*pow(c44,3) - 
51*c55*pow(c33,4)*pow(c44,3) + 4*c66*pow(c33,4)*pow(c44,3) + 15*pow(c33,5)*pow(c44,3) + 9*c33*c44*pow(c23,4)*pow(c55,2) + 
36*c22*c44*pow(c23,2)*pow(c33,2)*pow(c55,2) + 12*c44*c66*pow(c23,2)*pow(c33,2)*pow(c55,2) + 108*c44*pow(c23,3)*pow(c33,2)*pow(c55,2) + 
27*pow(c23,4)*pow(c33,2)*pow(c55,2) - 32*c12*c23*c44*pow(c33,3)*pow(c55,2) - 72*c22*c23*c44*pow(c33,3)*pow(c55,2) - 56*c23*c44*c66*pow(c33,3)*pow(c55,2) - 
36*c22*pow(c23,2)*pow(c33,3)*pow(c55,2) + 10*c44*pow(c23,2)*pow(c33,3)*pow(c55,2) - 12*c66*pow(c23,2)*pow(c33,3)*pow(c55,2) + 16*c12*c23*pow(c33,4)*pow(c55,2) + 
36*c11*c44*pow(c33,4)*pow(c55,2) + 16*c12*c44*pow(c33,4)*pow(c55,2) + 16*c23*c66*pow(c33,4)*pow(c55,2) + 28*c44*c66*pow(c33,4)*pow(c55,2) - 
12*c11*pow(c33,5)*pow(c55,2) + 17*c44*pow(c33,5)*pow(c55,2) - 4*c66*pow(c33,5)*pow(c55,2) - 3*pow(c33,6)*pow(c55,2) + 36*c33*pow(c23,3)*pow(c44,2)*pow(c55,2) + 
16*c12*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 72*c22*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 40*c23*c66*pow(c33,2)*pow(c44,2)*pow(c55,2) + 
170*pow(c23,2)*pow(c33,2)*pow(c44,2)*pow(c55,2) - 36*c11*pow(c33,3)*pow(c44,2)*pow(c55,2) - 32*c12*pow(c33,3)*pow(c44,2)*pow(c55,2) - 
36*c22*pow(c33,3)*pow(c44,2)*pow(c55,2) + 20*c23*pow(c33,3)*pow(c44,2)*pow(c55,2) - 56*c66*pow(c33,3)*pow(c44,2)*pow(c55,2) - 
42*pow(c33,4)*pow(c44,2)*pow(c55,2) + 36*c33*pow(c23,2)*pow(c44,3)*pow(c55,2) + 12*c11*pow(c33,2)*pow(c44,3)*pow(c55,2) + 
16*c12*pow(c33,2)*pow(c44,3)*pow(c55,2) + 36*c22*pow(c33,2)*pow(c44,3)*pow(c55,2) + 124*c23*pow(c33,2)*pow(c44,3)*pow(c55,2) + 
32*c66*pow(c33,2)*pow(c44,3)*pow(c55,2) + 64*pow(c33,3)*pow(c44,3)*pow(c55,2) - 
4*c13*(c33 - c44)*(-((c33 - c55)*c55*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2))) + 2*c12*c33*(c33 - c44)*(c23 + c44)*pow(c33 - c55,2) + 
2*c23*(c33 - c55)*(c33*c44*c55*(c44 + c55 + c66) - (c55*c66 + c44*(3*c55 + c66))*pow(c33,2) + c66*pow(c33,3) + pow(c44,2)*pow(c55,2)) + 
c33*(c33*c44*c55*(c44*(c55 + 6*c66) + c55*(19*c55 + 6*c66)) + 6*c11*(c33 - c55)*c55*pow(c33 - c44,2) + 
pow(c33,3)*(-(c44*c55) + 2*c44*c66 + 2*c55*c66 - 3*pow(c55,2)) - 4*(3*c55 + c66)*pow(c44,2)*pow(c55,2) - 
pow(c33,2)*((c55 + 2*c66)*pow(c44,2) + c44*(8*c55*c66 - 6*pow(c55,2)) + (9*c55 + 2*c66)*pow(c55,2)))) - 12*c22*c33*c44*pow(c23,2)*pow(c55,3) - 
4*c33*c44*c66*pow(c23,2)*pow(c55,3) - 36*c33*c44*pow(c23,3)*pow(c55,3) - 9*c33*pow(c23,4)*pow(c55,3) - 3*c44*pow(c23,4)*pow(c55,3) + 
16*c12*c23*c44*pow(c33,2)*pow(c55,3) + 24*c22*c23*c44*pow(c33,2)*pow(c55,3) + 24*c23*c44*c66*pow(c33,2)*pow(c55,3) + 12*c22*pow(c23,2)*pow(c33,2)*pow(c55,3) + 2*c44*pow(c23,2)*pow(c33,2)*pow(c55,3) + 4*c66*pow(c23,2)*pow(c33,2)*pow(c55,3) - 8*c12*c23*pow(c33,3)*pow(c55,3) - 36*c11*c44*pow(c33,3)*pow(c55,3) - 
8*c12*c44*pow(c33,3)*pow(c55,3) - 4*c23*c44*pow(c33,3)*pow(c55,3) - 8*c23*c66*pow(c33,3)*pow(c55,3) - 20*c44*c66*pow(c33,3)*pow(c55,3) - 
2*pow(c23,2)*pow(c33,3)*pow(c55,3) + 12*c11*pow(c33,4)*pow(c55,3) - 51*c44*pow(c33,4)*pow(c55,3) + 4*c66*pow(c33,4)*pow(c55,3) + 15*pow(c33,5)*pow(c55,3) - 
8*c12*c23*c33*pow(c44,2)*pow(c55,3) - 24*c22*c23*c33*pow(c44,2)*pow(c55,3) - 16*c23*c33*c66*pow(c44,2)*pow(c55,3) - 60*c33*pow(c23,2)*pow(c44,2)*pow(c55,3) - 12*pow(c23,3)*pow(c44,2)*pow(c55,3) + 36*c11*pow(c33,2)*pow(c44,2)*pow(c55,3) + 16*c12*pow(c33,2)*pow(c44,2)*pow(c55,3) + 
12*c22*pow(c33,2)*pow(c44,2)*pow(c55,3) + 4*c23*pow(c33,2)*pow(c44,2)*pow(c55,3) + 32*c66*pow(c33,2)*pow(c44,2)*pow(c55,3) + 64*pow(c33,3)*pow(c44,2)*pow(c55,3) - 
12*c11*c33*pow(c44,3)*pow(c55,3) - 8*c12*c33*pow(c44,3)*pow(c55,3) - 12*c22*c33*pow(c44,3)*pow(c55,3) - 48*c23*c33*pow(c44,3)*pow(c55,3) - 
16*c33*c66*pow(c44,3)*pow(c55,3) - 12*pow(c23,2)*pow(c44,3)*pow(c55,3) - 40*pow(c33,2)*pow(c44,3)*pow(c55,3) + 
2*(c33 - c44)*pow(c13,2)*(-2*c23*c44*(c33 - c55)*(c44*c55 + c33*(c44 + c55) - 3*pow(c33,2)) - 4*c44*c55*c66*pow(c33,2) + 
(c33 - c55)*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2)) - 6*c44*c55*pow(c33,3) + 4*c44*c66*pow(c33,3) + 2*c55*c66*pow(c33,3) + c44*pow(c33,4) + 3*c55*pow(c33,4) - 2*c66*pow(c33,4) - 6*c11*c33*(c33 - c55)*pow(c33 - c44,2) + 2*c33*c55*c66*pow(c44,2) - c55*pow(c33,2)*pow(c44,2) - 
2*c66*pow(c33,2)*pow(c44,2) + pow(c33,3)*pow(c44,2) - 55*c44*pow(c33,2)*pow(c55,2) + 27*pow(c33,3)*pow(c55,2) + 30*c33*pow(c44,2)*pow(c55,2) - 
12*c33*c44*pow(c55,3) + 6*pow(c33,2)*pow(c55,3) + 6*pow(c44,2)*pow(c55,3)) - 
3*cos(4*phi)*(4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + (3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) + 4*c44*(3*c33 + c44)*pow(c23,3)*pow(c33 - c55,3) + (3*c33 + c44)*pow(c23,4)*pow(c33 - c55,3) - 2*(c33 - c55)*pow(c23,2)*
((-c44 + c55 - 2*c66)*pow(c33,4) + pow(c33,3)*(2*c44*(c55 + c66) + c55*(c55 + 4*c66) - 9*pow(c44,2)) + 2*c33*c44*c55*(-3*c44*c55 + c55*c66 + 2*pow(c44,2)) + 2*c22*c33*(c33 - c44)*pow(c33 - c55,2) - 2*pow(c44,3)*pow(c55,2) - pow(c33,2)*(c44*c55*(5*c55 + 4*c66) - 17*c55*pow(c44,2) + 2*pow(c44,3) + 2*c66*pow(c55,2))
) - 4*c23*c33*(c33 - c55)*(-(c33*c44*c55*(5*c44*c55 + 2*c12*(c44 + 2*c55) + 6*c44*c66 + 6*c55*c66 - 5*pow(c44,2))) - 
pow(c33,3)*(-(c44*(c55 - 2*c66)) + 2*c55*(c12 + c66) + pow(c44,2)) + 2*c22*(c33 - c44)*c44*pow(c33 - c55,2) + 2*(c12 + 2*c66)*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(c44*c55*(4*c12 + c55 + 8*c66) + 2*(c55 + c66)*pow(c44,2) - 3*pow(c44,3) + 2*(c12 + c66)*pow(c55,2))) + 
4*c13*(c33 - c44)*(-((c33 - c55)*c55*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2))) + 2*c12*c33*(c33 - c44)*(c23 + c44)*pow(c33 - c55,2) + 
2*c23*(c33 - c55)*(c33*c44*c55*(c44 + c55 + c66) - (c55*c66 + c44*(3*c55 + c66))*pow(c33,2) + c66*pow(c33,3) + pow(c44,2)*pow(c55,2)) - 
c33*((c44*(c55 - 2*c66) - c55*(c55 + 2*c66))*pow(c33,3) + 2*c11*(c33 - c55)*c55*pow(c33 - c44,2) - 
c33*c44*c55*(5*c44*c55 + 6*c44*c66 + 6*c55*c66 - 5*pow(c55,2)) + 4*c66*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(2*c44*c55*(c55 + 4*c66) + (c55 + 2*c66)*pow(c44,2) + (-3*c55 + 2*c66)*pow(c55,2)))) - 
2*(c33 - c44)*pow(c13,2)*(-2*c23*c44*(c33 - c55)*(c44*c55 + c33*(c44 + c55) - 3*pow(c33,2)) - 4*c44*c55*c66*pow(c33,2) + 
(c33 - c55)*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2)) + 2*c44*c55*pow(c33,3) + 4*c44*c66*pow(c33,3) + 2*c55*c66*pow(c33,3) + c44*pow(c33,4) - c55*pow(c33,4) - 2*c66*pow(c33,4) + 2*c11*c33*(c33 - c55)*pow(c33 - c44,2) + 2*c33*c55*c66*pow(c44,2) - 5*c55*pow(c33,2)*pow(c44,2) - 
2*c66*pow(c33,2)*pow(c44,2) + pow(c33,3)*pow(c44,2) + 17*c44*pow(c33,2)*pow(c55,2) - 9*pow(c33,3)*pow(c55,2) - 6*c33*pow(c44,2)*pow(c55,2) + 
4*c33*c44*pow(c55,3) - 2*pow(c33,2)*pow(c55,3) - 2*pow(c44,2)*pow(c55,3)) - 
c33*(4*c22*(c33 - c44)*pow(c44,2)*pow(c33 - c55,3) + pow(c33,5)*pow(c44 - c55,2) - 
4*c33*(c11*(c44 + 3*c55) - 4*(c44 + c55)*(c12 + 2*c66))*pow(c44,2)*pow(c55,2) - 
pow(c33,4)*(c44*c55*(8*c12 - 5*c55 + 8*c66) + (-5*c55 + 4*c66)*pow(c44,2) + 5*pow(c44,3) + (-4*c11 + 5*c55 + 4*c66)*pow(c55,2)) - 
4*c44*c55*pow(c33,2)*(-3*c11*c55*(c44 + c55) + 2*c12*(4*c44*c55 + pow(c44,2) + pow(c55,2)) + c66*(14*c44*c55 + 5*pow(c44,2) + 5*pow(c55,2))) + 
4*(c11 - 2*(c12 + 2*c66))*pow(c44,3)*pow(c55,3) + pow(c33,3)*
(16*c12*c44*c55*(c44 + c55) + 2*c55*(-9*c55 + 14*c66)*pow(c44,2) + (9*c55 + 4*c66)*pow(c44,3) + c44*(-12*c11 + 9*c55 + 28*c66)*pow(c55,2) + 
4*(-c11 + c66)*pow(c55,3))))))/8.;
	
	
	
	ppsi13[j] = (-3*pow(c33,-0.5)*pow(c33 - c44,-3)*pow(c33 - c55,-3)*(-3*c44*c55*pow(c23,4)*pow(c33,2) - 12*c22*c44*c55*pow(c23,2)*pow(c33,3) - 36*c44*c55*pow(c23,3)*pow(c33,3) + 
c44*pow(c23,4)*pow(c33,3) - 9*c55*pow(c23,4)*pow(c33,3) + 24*c22*c23*c44*c55*pow(c33,4) + 4*c22*c44*pow(c23,2)*pow(c33,4) + 12*c22*c55*pow(c23,2)*pow(c33,4) - 6*c44*c55*pow(c23,2)*pow(c33,4) + 12*c44*pow(c23,3)*pow(c33,4) + 3*pow(c23,4)*pow(c33,4) - 8*c22*c23*c44*pow(c33,5) - 4*c22*pow(c23,2)*pow(c33,5) + 
2*c44*pow(c23,2)*pow(c33,5) - 4*c13*c33*c55*(-2*c11*(c33 - c55) + c55*(c33 + 3*c55))*pow(c33 - c44,3) - 4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) - 
(3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) - 12*c55*pow(c23,3)*pow(c33,2)*pow(c44,2) - 24*c22*c23*c55*pow(c33,3)*pow(c44,2) - 
54*c55*pow(c23,2)*pow(c33,3)*pow(c44,2) + 4*pow(c23,3)*pow(c33,3)*pow(c44,2) + 8*c22*c23*pow(c33,4)*pow(c44,2) + 12*c22*c55*pow(c33,4)*pow(c44,2) - 
12*c23*c55*pow(c33,4)*pow(c44,2) + 18*pow(c23,2)*pow(c33,4)*pow(c44,2) - 4*c22*pow(c33,5)*pow(c44,2) + 4*c23*pow(c33,5)*pow(c44,2) + 3*c55*pow(c33,5)*pow(c44,2) - 
pow(c33,6)*pow(c44,2) - 12*c55*pow(c23,2)*pow(c33,2)*pow(c44,3) - 12*c22*c55*pow(c33,3)*pow(c44,3) - 36*c23*c55*pow(c33,3)*pow(c44,3) + 
4*pow(c23,2)*pow(c33,3)*pow(c44,3) + 4*c22*pow(c33,4)*pow(c44,3) + 12*c23*pow(c33,4)*pow(c44,3) - 15*c55*pow(c33,4)*pow(c44,3) + 5*pow(c33,5)*pow(c44,3) + 
3*c33*c44*pow(c23,4)*pow(c55,2) + 12*c22*c44*pow(c23,2)*pow(c33,2)*pow(c55,2) + 36*c44*pow(c23,3)*pow(c33,2)*pow(c55,2) + 9*pow(c23,4)*pow(c33,2)*pow(c55,2) - 24*c22*c23*c44*pow(c33,3)*pow(c55,2) - 12*c22*pow(c23,2)*pow(c33,3)*pow(c55,2) + 6*c44*pow(c23,2)*pow(c33,3)*pow(c55,2) - 12*c11*c44*pow(c33,4)*pow(c55,2) + 4*c11*pow(c33,5)*pow(c55,2) - 3*c44*pow(c33,5)*pow(c55,2) + pow(c33,6)*pow(c55,2) + 12*c33*pow(c23,3)*pow(c44,2)*pow(c55,2) + 
24*c22*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 54*pow(c23,2)*pow(c33,2)*pow(c44,2)*pow(c55,2) + 12*c11*pow(c33,3)*pow(c44,2)*pow(c55,2) - 
12*c22*pow(c33,3)*pow(c44,2)*pow(c55,2) + 12*c23*pow(c33,3)*pow(c44,2)*pow(c55,2) + 12*c33*pow(c23,2)*pow(c44,3)*pow(c55,2) - 
4*c11*pow(c33,2)*pow(c44,3)*pow(c55,2) + 12*c22*pow(c33,2)*pow(c44,3)*pow(c55,2) + 36*c23*pow(c33,2)*pow(c44,3)*pow(c55,2) + 14*pow(c33,3)*pow(c44,3)*pow(c55,2) + 
2*pow(c13,2)*pow(c33 - c44,3)*(2*c11*c33*(c33 - c55) - c55*(9*c33*c55 + pow(c33,2) + 2*pow(c55,2))) - 4*c22*c33*c44*pow(c23,2)*pow(c55,3) - 
12*c33*c44*pow(c23,3)*pow(c55,3) - 3*c33*pow(c23,4)*pow(c55,3) - c44*pow(c23,4)*pow(c55,3) + 8*c22*c23*c44*pow(c33,2)*pow(c55,3) + 
4*c22*pow(c23,2)*pow(c33,2)*pow(c55,3) - 2*c44*pow(c23,2)*pow(c33,2)*pow(c55,3) + 12*c11*c44*pow(c33,3)*pow(c55,3) - 4*c11*pow(c33,4)*pow(c55,3) + 
15*c44*pow(c33,4)*pow(c55,3) - 5*pow(c33,5)*pow(c55,3) - 8*c22*c23*c33*pow(c44,2)*pow(c55,3) - 18*c33*pow(c23,2)*pow(c44,2)*pow(c55,3) - 
4*pow(c23,3)*pow(c44,2)*pow(c55,3) - 12*c11*pow(c33,2)*pow(c44,2)*pow(c55,3) + 4*c22*pow(c33,2)*pow(c44,2)*pow(c55,3) - 4*c23*pow(c33,2)*pow(c44,2)*pow(c55,3) - 
14*pow(c33,3)*pow(c44,2)*pow(c55,3) + 4*c11*c33*pow(c44,3)*pow(c55,3) - 4*c22*c33*pow(c44,3)*pow(c55,3) - 12*c23*c33*pow(c44,3)*pow(c55,3) - 
4*pow(c23,2)*pow(c44,3)*pow(c55,3) + cos(2*phi)*(4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + (3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) + 
4*c44*(3*c33 + c44)*pow(c23,3)*pow(c33 - c55,3) + (3*c33 + c44)*pow(c23,4)*pow(c33 - c55,3) - 
2*(c33 - c55)*pow(c23,2)*((-c44 + c55 - 2*c66)*pow(c33,4) + pow(c33,3)*(2*c44*(c55 + c66) + c55*(c55 + 4*c66) - 9*pow(c44,2)) + 
2*c33*c44*c55*(-3*c44*c55 + c55*c66 + 2*pow(c44,2)) + 2*c22*c33*(c33 - c44)*pow(c33 - c55,2) - 2*pow(c44,3)*pow(c55,2) - 
pow(c33,2)*(c44*c55*(5*c55 + 4*c66) - 17*c55*pow(c44,2) + 2*pow(c44,3) + 2*c66*pow(c55,2))) - 
4*c23*c33*(c33 - c55)*(-(c33*c44*c55*(5*c44*c55 + 2*c12*(c44 + 2*c55) + 6*c44*c66 + 6*c55*c66 - 5*pow(c44,2))) - 
pow(c33,3)*(-(c44*(c55 - 2*c66)) + 2*c55*(c12 + c66) + pow(c44,2)) + 2*c22*(c33 - c44)*c44*pow(c33 - c55,2) + 2*(c12 + 2*c66)*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(c44*c55*(4*c12 + c55 + 8*c66) + 2*(c55 + c66)*pow(c44,2) - 3*pow(c44,3) + 2*(c12 + c66)*pow(c55,2))) + 
4*c13*(c33 - c44)*(-((c33 - c55)*c55*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2))) + 2*c12*c33*(c33 - c44)*(c23 + c44)*pow(c33 - c55,2) + 
2*c23*(c33 - c55)*(c33*c44*c55*(c44 + c55 + c66) - (c55*c66 + c44*(3*c55 + c66))*pow(c33,2) + c66*pow(c33,3) + pow(c44,2)*pow(c55,2)) - 
c33*((c44*(c55 - 2*c66) - c55*(c55 + 2*c66))*pow(c33,3) + 2*c11*(c33 - c55)*c55*pow(c33 - c44,2) - 
c33*c44*c55*(5*c44*c55 + 6*c44*c66 + 6*c55*c66 - 5*pow(c55,2)) + 4*c66*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(2*c44*c55*(c55 + 4*c66) + (c55 + 2*c66)*pow(c44,2) + (-3*c55 + 2*c66)*pow(c55,2)))) - 
2*(c33 - c44)*pow(c13,2)*(-2*c23*c44*(c33 - c55)*(c44*c55 + c33*(c44 + c55) - 3*pow(c33,2)) - 4*c44*c55*c66*pow(c33,2) + 
(c33 - c55)*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2)) + 2*c44*c55*pow(c33,3) + 4*c44*c66*pow(c33,3) + 2*c55*c66*pow(c33,3) + c44*pow(c33,4) - c55*pow(c33,4) - 2*c66*pow(c33,4) + 2*c11*c33*(c33 - c55)*pow(c33 - c44,2) + 2*c33*c55*c66*pow(c44,2) - 5*c55*pow(c33,2)*pow(c44,2) - 
2*c66*pow(c33,2)*pow(c44,2) + pow(c33,3)*pow(c44,2) + 17*c44*pow(c33,2)*pow(c55,2) - 9*pow(c33,3)*pow(c55,2) - 6*c33*pow(c44,2)*pow(c55,2) + 
4*c33*c44*pow(c55,3) - 2*pow(c33,2)*pow(c55,3) - 2*pow(c44,2)*pow(c55,3)) - 
c33*(4*c22*(c33 - c44)*pow(c44,2)*pow(c33 - c55,3) + pow(c33,5)*pow(c44 - c55,2) - 
4*c33*(c11*(c44 + 3*c55) - 4*(c44 + c55)*(c12 + 2*c66))*pow(c44,2)*pow(c55,2) - 
pow(c33,4)*(c44*c55*(8*c12 - 5*c55 + 8*c66) + (-5*c55 + 4*c66)*pow(c44,2) + 5*pow(c44,3) + (-4*c11 + 5*c55 + 4*c66)*pow(c55,2)) - 
4*c44*c55*pow(c33,2)*(-3*c11*c55*(c44 + c55) + 2*c12*(4*c44*c55 + pow(c44,2) + pow(c55,2)) + c66*(14*c44*c55 + 5*pow(c44,2) + 5*pow(c55,2))) + 
4*(c11 - 2*(c12 + 2*c66))*pow(c44,3)*pow(c55,3) + pow(c33,3)*
(16*c12*c44*c55*(c44 + c55) + 2*c55*(-9*c55 + 14*c66)*pow(c44,2) + (9*c55 + 4*c66)*pow(c44,3) + c44*(-12*c11 + 9*c55 + 28*c66)*pow(c55,2) + 
4*(-c11 + c66)*pow(c55,3)))))*sin(2*phi))/4.;



	ppsi04[j] = (3*pow(c33,-0.5)*pow(c33 - c44,-3)*pow(c33 - c55,-3)*(-9*c44*c55*pow(c13,4)*pow(c33,2) + 8*c44*c55*pow(c13,2)*pow(c23,2)*pow(c33,2) - 9*c44*c55*pow(c23,4)*pow(c33,2) - 
32*c12*c13*c23*c44*c55*pow(c33,3) - 32*c13*c23*c44*c55*c66*pow(c33,3) - 36*c11*c44*c55*pow(c13,2)*pow(c33,3) - 16*c23*c44*c55*pow(c13,2)*pow(c33,3) - 
12*c44*c55*c66*pow(c13,2)*pow(c33,3) - 108*c44*c55*pow(c13,3)*pow(c33,3) - 27*c44*pow(c13,4)*pow(c33,3) + 3*c55*pow(c13,4)*pow(c33,3) - 
16*c13*c44*c55*pow(c23,2)*pow(c33,3) - 36*c22*c44*c55*pow(c23,2)*pow(c33,3) - 12*c44*c55*c66*pow(c23,2)*pow(c33,3) - 8*c44*pow(c13,2)*pow(c23,2)*pow(c33,3) - 
8*c55*pow(c13,2)*pow(c23,2)*pow(c33,3) - 108*c44*c55*pow(c23,3)*pow(c33,3) + 3*c44*pow(c23,4)*pow(c33,3) - 27*c55*pow(c23,4)*pow(c33,3) + 
16*c12*c13*c23*c44*pow(c33,4) + 16*c12*c13*c23*c55*pow(c33,4) + 72*c11*c13*c44*c55*pow(c33,4) + 16*c12*c13*c44*c55*pow(c33,4) + 16*c12*c23*c44*c55*pow(c33,4) + 
24*c13*c23*c44*c55*pow(c33,4) + 72*c22*c23*c44*c55*pow(c33,4) + 16*c13*c23*c44*c66*pow(c33,4) + 16*c13*c23*c55*c66*pow(c33,4) + 40*c13*c44*c55*c66*pow(c33,4) + 
40*c23*c44*c55*c66*pow(c33,4) + 36*c11*c44*pow(c13,2)*pow(c33,4) + 12*c23*c44*pow(c13,2)*pow(c33,4) + 12*c11*c55*pow(c13,2)*pow(c33,4) - 
18*c44*c55*pow(c13,2)*pow(c33,4) + 12*c44*c66*pow(c13,2)*pow(c33,4) + 4*c55*c66*pow(c13,2)*pow(c33,4) + 36*c55*pow(c13,3)*pow(c33,4) + 9*pow(c13,4)*pow(c33,4) + 
12*c22*c44*pow(c23,2)*pow(c33,4) + 12*c13*c55*pow(c23,2)*pow(c33,4) + 36*c22*c55*pow(c23,2)*pow(c33,4) - 18*c44*c55*pow(c23,2)*pow(c33,4) + 
4*c44*c66*pow(c23,2)*pow(c33,4) + 12*c55*c66*pow(c23,2)*pow(c33,4) + 6*pow(c13,2)*pow(c23,2)*pow(c33,4) + 36*c44*pow(c23,3)*pow(c33,4) + 9*pow(c23,4)*pow(c33,4) - 
8*c12*c13*c23*pow(c33,5) - 8*c12*c13*c44*pow(c33,5) - 24*c22*c23*c44*pow(c33,5) - 24*c11*c13*c55*pow(c33,5) - 8*c12*c23*c55*pow(c33,5) - 
8*c12*c44*c55*pow(c33,5) + 4*c13*c44*c55*pow(c33,5) + 4*c23*c44*c55*pow(c33,5) - 8*c13*c23*c66*pow(c33,5) - 8*c13*c44*c66*pow(c33,5) - 8*c23*c44*c66*pow(c33,5) - 
8*c13*c55*c66*pow(c33,5) - 8*c23*c55*c66*pow(c33,5) - 8*c44*c55*c66*pow(c33,5) - 12*c11*pow(c13,2)*pow(c33,5) + 2*c44*pow(c13,2)*pow(c33,5) + 
6*c55*pow(c13,2)*pow(c33,5) - 4*c66*pow(c13,2)*pow(c33,5) - 12*c22*pow(c23,2)*pow(c33,5) + 6*c44*pow(c23,2)*pow(c33,5) + 2*c55*pow(c23,2)*pow(c33,5) - 
4*c66*pow(c23,2)*pow(c33,5) - 2*c44*c55*pow(c33,6) + 9*c33*c55*pow(c13,4)*pow(c44,2) + 16*c12*c13*c23*c55*pow(c33,2)*pow(c44,2) + 
16*c13*c23*c55*c66*pow(c33,2)*pow(c44,2) + 36*c11*c55*pow(c13,2)*pow(c33,2)*pow(c44,2) + 16*c23*c55*pow(c13,2)*pow(c33,2)*pow(c44,2) + 
12*c55*c66*pow(c13,2)*pow(c33,2)*pow(c44,2) + 108*c55*pow(c13,3)*pow(c33,2)*pow(c44,2) + 27*pow(c13,4)*pow(c33,2)*pow(c44,2) + 
4*c13*c55*pow(c23,2)*pow(c33,2)*pow(c44,2) + 2*pow(c13,2)*pow(c23,2)*pow(c33,2)*pow(c44,2) - 36*c55*pow(c23,3)*pow(c33,2)*pow(c44,2) - 
8*c12*c13*c23*pow(c33,3)*pow(c44,2) - 72*c11*c13*c55*pow(c33,3)*pow(c44,2) - 32*c12*c13*c55*pow(c33,3)*pow(c44,2) - 8*c12*c23*c55*pow(c33,3)*pow(c44,2) - 
32*c13*c23*c55*pow(c33,3)*pow(c44,2) - 72*c22*c23*c55*pow(c33,3)*pow(c44,2) - 8*c13*c23*c66*pow(c33,3)*pow(c44,2) - 56*c13*c55*c66*pow(c33,3)*pow(c44,2) - 
32*c23*c55*c66*pow(c33,3)*pow(c44,2) - 36*c11*pow(c13,2)*pow(c33,3)*pow(c44,2) - 16*c23*pow(c13,2)*pow(c33,3)*pow(c44,2) + 
10*c55*pow(c13,2)*pow(c33,3)*pow(c44,2) - 12*c66*pow(c13,2)*pow(c33,3)*pow(c44,2) - 164*c55*pow(c23,2)*pow(c33,3)*pow(c44,2) + 
12*pow(c23,3)*pow(c33,3)*pow(c44,2) + 16*c12*c13*pow(c33,4)*pow(c44,2) + 24*c22*c23*pow(c33,4)*pow(c44,2) + 16*c12*c55*pow(c33,4)*pow(c44,2) + 
36*c22*c55*pow(c33,4)*pow(c44,2) - 36*c23*c55*pow(c33,4)*pow(c44,2) + 16*c13*c66*pow(c33,4)*pow(c44,2) + 8*c23*c66*pow(c33,4)*pow(c44,2) + 
28*c55*c66*pow(c33,4)*pow(c44,2) + 54*pow(c23,2)*pow(c33,4)*pow(c44,2) - 12*c22*pow(c33,5)*pow(c44,2) + 12*c23*pow(c33,5)*pow(c44,2) + 
17*c55*pow(c33,5)*pow(c44,2) - 4*c66*pow(c33,5)*pow(c44,2) - 3*pow(c33,6)*pow(c44,2) - 12*c11*c33*c55*pow(c13,2)*pow(c44,3) - 
4*c33*c55*c66*pow(c13,2)*pow(c44,3) - 36*c33*c55*pow(c13,3)*pow(c44,3) - 9*c33*pow(c13,4)*pow(c44,3) - 3*c55*pow(c13,4)*pow(c44,3) + 
24*c11*c13*c55*pow(c33,2)*pow(c44,3) + 16*c12*c13*c55*pow(c33,2)*pow(c44,3) + 8*c13*c23*c55*pow(c33,2)*pow(c44,3) + 24*c13*c55*c66*pow(c33,2)*pow(c44,3) + 
12*c11*pow(c13,2)*pow(c33,2)*pow(c44,3) + 4*c23*pow(c13,2)*pow(c33,2)*pow(c44,3) + 2*c55*pow(c13,2)*pow(c33,2)*pow(c44,3) + 
4*c66*pow(c13,2)*pow(c33,2)*pow(c44,3) - 36*c55*pow(c23,2)*pow(c33,2)*pow(c44,3) - 8*c12*c13*pow(c33,3)*pow(c44,3) - 8*c12*c55*pow(c33,3)*pow(c44,3) - 
4*c13*c55*pow(c33,3)*pow(c44,3) - 36*c22*c55*pow(c33,3)*pow(c44,3) - 112*c23*c55*pow(c33,3)*pow(c44,3) - 8*c13*c66*pow(c33,3)*pow(c44,3) - 
20*c55*c66*pow(c33,3)*pow(c44,3) - 2*pow(c13,2)*pow(c33,3)*pow(c44,3) + 12*pow(c23,2)*pow(c33,3)*pow(c44,3) + 12*c22*pow(c33,4)*pow(c44,3) + 
36*c23*pow(c33,4)*pow(c44,3) - 51*c55*pow(c33,4)*pow(c44,3) + 4*c66*pow(c33,4)*pow(c44,3) + 15*pow(c33,5)*pow(c44,3) + 9*c33*c44*pow(c23,4)*pow(c55,2) + 
16*c12*c13*c23*c44*pow(c33,2)*pow(c55,2) + 16*c13*c23*c44*c66*pow(c33,2)*pow(c55,2) + 4*c23*c44*pow(c13,2)*pow(c33,2)*pow(c55,2) - 
36*c44*pow(c13,3)*pow(c33,2)*pow(c55,2) + 16*c13*c44*pow(c23,2)*pow(c33,2)*pow(c55,2) + 36*c22*c44*pow(c23,2)*pow(c33,2)*pow(c55,2) + 
12*c44*c66*pow(c23,2)*pow(c33,2)*pow(c55,2) + 2*pow(c13,2)*pow(c23,2)*pow(c33,2)*pow(c55,2) + 108*c44*pow(c23,3)*pow(c33,2)*pow(c55,2) + 
27*pow(c23,4)*pow(c33,2)*pow(c55,2) - 8*c12*c13*c23*pow(c33,3)*pow(c55,2) - 72*c11*c13*c44*pow(c33,3)*pow(c55,2) - 8*c12*c13*c44*pow(c33,3)*pow(c55,2) - 
32*c12*c23*c44*pow(c33,3)*pow(c55,2) - 32*c13*c23*c44*pow(c33,3)*pow(c55,2) - 72*c22*c23*c44*pow(c33,3)*pow(c55,2) - 8*c13*c23*c66*pow(c33,3)*pow(c55,2) - 
32*c13*c44*c66*pow(c33,3)*pow(c55,2) - 56*c23*c44*c66*pow(c33,3)*pow(c55,2) - 164*c44*pow(c13,2)*pow(c33,3)*pow(c55,2) + 12*pow(c13,3)*pow(c33,3)*pow(c55,2) - 16*c13*pow(c23,2)*pow(c33,3)*pow(c55,2) - 36*c22*pow(c23,2)*pow(c33,3)*pow(c55,2) + 10*c44*pow(c23,2)*pow(c33,3)*pow(c55,2) - 
12*c66*pow(c23,2)*pow(c33,3)*pow(c55,2) + 24*c11*c13*pow(c33,4)*pow(c55,2) + 16*c12*c23*pow(c33,4)*pow(c55,2) + 36*c11*c44*pow(c33,4)*pow(c55,2) + 
16*c12*c44*pow(c33,4)*pow(c55,2) - 36*c13*c44*pow(c33,4)*pow(c55,2) + 8*c13*c66*pow(c33,4)*pow(c55,2) + 16*c23*c66*pow(c33,4)*pow(c55,2) + 
28*c44*c66*pow(c33,4)*pow(c55,2) + 54*pow(c13,2)*pow(c33,4)*pow(c55,2) - 12*c11*pow(c33,5)*pow(c55,2) + 12*c13*pow(c33,5)*pow(c55,2) + 
17*c44*pow(c33,5)*pow(c55,2) - 4*c66*pow(c33,5)*pow(c55,2) - 3*pow(c33,6)*pow(c55,2) - 8*c12*c13*c23*c33*pow(c44,2)*pow(c55,2) - 
8*c13*c23*c33*c66*pow(c44,2)*pow(c55,2) + 36*c33*pow(c13,3)*pow(c44,2)*pow(c55,2) - 2*pow(c13,2)*pow(c23,2)*pow(c44,2)*pow(c55,2) + 
36*c33*pow(c23,3)*pow(c44,2)*pow(c55,2) + 72*c11*c13*pow(c33,2)*pow(c44,2)*pow(c55,2) + 16*c12*c13*pow(c33,2)*pow(c44,2)*pow(c55,2) + 
16*c12*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 32*c13*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 72*c22*c23*pow(c33,2)*pow(c44,2)*pow(c55,2) + 
40*c13*c66*pow(c33,2)*pow(c44,2)*pow(c55,2) + 40*c23*c66*pow(c33,2)*pow(c44,2)*pow(c55,2) + 170*pow(c13,2)*pow(c33,2)*pow(c44,2)*pow(c55,2) + 
170*pow(c23,2)*pow(c33,2)*pow(c44,2)*pow(c55,2) - 36*c11*pow(c33,3)*pow(c44,2)*pow(c55,2) - 32*c12*pow(c33,3)*pow(c44,2)*pow(c55,2) + 
20*c13*pow(c33,3)*pow(c44,2)*pow(c55,2) - 36*c22*pow(c33,3)*pow(c44,2)*pow(c55,2) + 20*c23*pow(c33,3)*pow(c44,2)*pow(c55,2) - 
56*c66*pow(c33,3)*pow(c44,2)*pow(c55,2) - 42*pow(c33,4)*pow(c44,2)*pow(c55,2) - 24*c11*c13*c33*pow(c44,3)*pow(c55,2) - 8*c12*c13*c33*pow(c44,3)*pow(c55,2) - 16*c13*c33*c66*pow(c44,3)*pow(c55,2) - 4*c23*pow(c13,2)*pow(c44,3)*pow(c55,2) - 60*c33*pow(c13,2)*pow(c44,3)*pow(c55,2) - 12*pow(c13,3)*pow(c44,3)*pow(c55,2) + 36*c33*pow(c23,2)*pow(c44,3)*pow(c55,2) + 12*c11*pow(c33,2)*pow(c44,3)*pow(c55,2) + 16*c12*pow(c33,2)*pow(c44,3)*pow(c55,2) + 
4*c13*pow(c33,2)*pow(c44,3)*pow(c55,2) + 36*c22*pow(c33,2)*pow(c44,3)*pow(c55,2) + 124*c23*pow(c33,2)*pow(c44,3)*pow(c55,2) + 
32*c66*pow(c33,2)*pow(c44,3)*pow(c55,2) + 64*pow(c33,3)*pow(c44,3)*pow(c55,2) - 12*c22*c33*c44*pow(c23,2)*pow(c55,3) - 4*c33*c44*c66*pow(c23,2)*pow(c55,3) - 36*c33*c44*pow(c23,3)*pow(c55,3) - 9*c33*pow(c23,4)*pow(c55,3) - 3*c44*pow(c23,4)*pow(c55,3) + 16*c12*c23*c44*pow(c33,2)*pow(c55,3) + 
8*c13*c23*c44*pow(c33,2)*pow(c55,3) + 24*c22*c23*c44*pow(c33,2)*pow(c55,3) + 24*c23*c44*c66*pow(c33,2)*pow(c55,3) - 36*c44*pow(c13,2)*pow(c33,2)*pow(c55,3) + 4*c13*pow(c23,2)*pow(c33,2)*pow(c55,3) + 12*c22*pow(c23,2)*pow(c33,2)*pow(c55,3) + 2*c44*pow(c23,2)*pow(c33,2)*pow(c55,3) + 
4*c66*pow(c23,2)*pow(c33,2)*pow(c55,3) - 8*c12*c23*pow(c33,3)*pow(c55,3) - 36*c11*c44*pow(c33,3)*pow(c55,3) - 8*c12*c44*pow(c33,3)*pow(c55,3) - 
112*c13*c44*pow(c33,3)*pow(c55,3) - 4*c23*c44*pow(c33,3)*pow(c55,3) - 8*c23*c66*pow(c33,3)*pow(c55,3) - 20*c44*c66*pow(c33,3)*pow(c55,3) + 
12*pow(c13,2)*pow(c33,3)*pow(c55,3) - 2*pow(c23,2)*pow(c33,3)*pow(c55,3) + 12*c11*pow(c33,4)*pow(c55,3) + 36*c13*pow(c33,4)*pow(c55,3) - 
51*c44*pow(c33,4)*pow(c55,3) + 4*c66*pow(c33,4)*pow(c55,3) + 15*pow(c33,5)*pow(c55,3) - 8*c12*c23*c33*pow(c44,2)*pow(c55,3) - 
24*c22*c23*c33*pow(c44,2)*pow(c55,3) - 16*c23*c33*c66*pow(c44,2)*pow(c55,3) + 36*c33*pow(c13,2)*pow(c44,2)*pow(c55,3) - 4*c13*pow(c23,2)*pow(c44,2)*pow(c55,3) - 60*c33*pow(c23,2)*pow(c44,2)*pow(c55,3) - 12*pow(c23,3)*pow(c44,2)*pow(c55,3) + 36*c11*pow(c33,2)*pow(c44,2)*pow(c55,3) + 
16*c12*pow(c33,2)*pow(c44,2)*pow(c55,3) + 124*c13*pow(c33,2)*pow(c44,2)*pow(c55,3) + 12*c22*pow(c33,2)*pow(c44,2)*pow(c55,3) + 
4*c23*pow(c33,2)*pow(c44,2)*pow(c55,3) + 32*c66*pow(c33,2)*pow(c44,2)*pow(c55,3) + 64*pow(c33,3)*pow(c44,2)*pow(c55,3) - 8*c13*c23*pow(c44,3)*pow(c55,3) - 
12*c11*c33*pow(c44,3)*pow(c55,3) - 8*c12*c33*pow(c44,3)*pow(c55,3) - 48*c13*c33*pow(c44,3)*pow(c55,3) - 12*c22*c33*pow(c44,3)*pow(c55,3) - 
48*c23*c33*pow(c44,3)*pow(c55,3) - 16*c33*c66*pow(c44,3)*pow(c55,3) - 12*pow(c13,2)*pow(c44,3)*pow(c55,3) - 12*pow(c23,2)*pow(c44,3)*pow(c55,3) - 
40*pow(c33,2)*pow(c44,3)*pow(c55,3) - 4*cos(2*phi)*(4*c13*c33*c55*(-2*c11*(c33 - c55) + c55*(c33 + 3*c55))*pow(c33 - c44,3) + 
4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + (3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) - 
4*c23*c33*c44*(-2*c22*(c33 - c44) + c44*(c33 + 3*c44))*pow(c33 - c55,3) - 4*c44*(3*c33 + c44)*pow(c23,3)*pow(c33 - c55,3) - 
(3*c33 + c44)*pow(c23,4)*pow(c33 - c55,3) + 2*pow(c23,2)*(2*c22*c33*(c33 - c44) - c44*(9*c33*c44 + pow(c33,2) + 2*pow(c44,2)))*pow(c33 - c55,3) - 
2*pow(c13,2)*pow(c33 - c44,3)*(2*c11*c33*(c33 - c55) - c55*(9*c33*c55 + pow(c33,2) + 2*pow(c55,2))) + 
c33*(4*c22*(c33 - c44)*pow(c44,2)*pow(c33 - c55,3) + pow(c33,5)*(pow(c44,2) - pow(c55,2)) - 
2*c44*(7*c44*(c44 - c55) + 6*c11*(c44 + c55))*pow(c33,2)*pow(c55,2) + 4*c11*c33*(c44 + 3*c55)*pow(c44,2)*pow(c55,2) + 
c55*pow(c33,3)*(3*c44*(4*c11 - 5*c55)*c55 + 15*pow(c44,3) + 4*c11*pow(c55,2)) + 
pow(c33,4)*(-3*c55*pow(c44,2) - 5*pow(c44,3) + 3*c44*pow(c55,2) + (-4*c11 + 5*c55)*pow(c55,2)) - 4*c11*pow(c44,3)*pow(c55,3))) + 
cos(4*phi)*(4*c55*(3*c33 + c55)*pow(c13,3)*pow(c33 - c44,3) + (3*c33 + c55)*pow(c13,4)*pow(c33 - c44,3) + 4*c44*(3*c33 + c44)*pow(c23,3)*pow(c33 - c55,3) + 
(3*c33 + c44)*pow(c23,4)*pow(c33 - c55,3) - 2*(c33 - c55)*pow(c23,2)*
((-c44 + c55 - 2*c66)*pow(c33,4) + pow(c33,3)*(2*c44*(c55 + c66) + c55*(c55 + 4*c66) - 9*pow(c44,2)) + 2*c33*c44*c55*(-3*c44*c55 + c55*c66 + 2*pow(c44,2)) + 2*c22*c33*(c33 - c44)*pow(c33 - c55,2) - 2*pow(c44,3)*pow(c55,2) - pow(c33,2)*(c44*c55*(5*c55 + 4*c66) - 17*c55*pow(c44,2) + 2*pow(c44,3) + 2*c66*pow(c55,2))) - 4*c23*c33*(c33 - c55)*(-(c33*c44*c55*(5*c44*c55 + 2*c12*(c44 + 2*c55) + 6*c44*c66 + 6*c55*c66 - 5*pow(c44,2))) - 
pow(c33,3)*(-(c44*(c55 - 2*c66)) + 2*c55*(c12 + c66) + pow(c44,2)) + 2*c22*(c33 - c44)*c44*pow(c33 - c55,2) + 2*(c12 + 2*c66)*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(c44*c55*(4*c12 + c55 + 8*c66) + 2*(c55 + c66)*pow(c44,2) - 3*pow(c44,3) + 2*(c12 + c66)*pow(c55,2))) + 
4*c13*(c33 - c44)*(-((c33 - c55)*c55*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2))) + 2*c12*c33*(c33 - c44)*(c23 + c44)*pow(c33 - c55,2) + 
2*c23*(c33 - c55)*(c33*c44*c55*(c44 + c55 + c66) - (c55*c66 + c44*(3*c55 + c66))*pow(c33,2) + c66*pow(c33,3) + pow(c44,2)*pow(c55,2)) - 
c33*((c44*(c55 - 2*c66) - c55*(c55 + 2*c66))*pow(c33,3) + 2*c11*(c33 - c55)*c55*pow(c33 - c44,2) - 
c33*c44*c55*(5*c44*c55 + 6*c44*c66 + 6*c55*c66 - 5*pow(c55,2)) + 4*c66*pow(c44,2)*pow(c55,2) + 
pow(c33,2)*(2*c44*c55*(c55 + 4*c66) + (c55 + 2*c66)*pow(c44,2) + (-3*c55 + 2*c66)*pow(c55,2)))) - 
2*(c33 - c44)*pow(c13,2)*(-2*c23*c44*(c33 - c55)*(c44*c55 + c33*(c44 + c55) - 3*pow(c33,2)) - 4*c44*c55*c66*pow(c33,2) + 
(c33 - c55)*pow(c23,2)*(-(c44*c55) - c33*(c44 + c55) + 3*pow(c33,2)) + 2*c44*c55*pow(c33,3) + 4*c44*c66*pow(c33,3) + 2*c55*c66*pow(c33,3) + c44*pow(c33,4) - c55*pow(c33,4) - 2*c66*pow(c33,4) + 2*c11*c33*(c33 - c55)*pow(c33 - c44,2) + 2*c33*c55*c66*pow(c44,2) - 5*c55*pow(c33,2)*pow(c44,2) - 
2*c66*pow(c33,2)*pow(c44,2) + pow(c33,3)*pow(c44,2) + 17*c44*pow(c33,2)*pow(c55,2) - 9*pow(c33,3)*pow(c55,2) - 6*c33*pow(c44,2)*pow(c55,2) + 
4*c33*c44*pow(c55,3) - 2*pow(c33,2)*pow(c55,3) - 2*pow(c44,2)*pow(c55,3)) - 
c33*(4*c22*(c33 - c44)*pow(c44,2)*pow(c33 - c55,3) + pow(c33,5)*pow(c44 - c55,2) - 
4*c33*(c11*(c44 + 3*c55) - 4*(c44 + c55)*(c12 + 2*c66))*pow(c44,2)*pow(c55,2) - 
pow(c33,4)*(c44*c55*(8*c12 - 5*c55 + 8*c66) + (-5*c55 + 4*c66)*pow(c44,2) + 5*pow(c44,3) + (-4*c11 + 5*c55 + 4*c66)*pow(c55,2)) - 
4*c44*c55*pow(c33,2)*(-3*c11*c55*(c44 + c55) + 2*c12*(4*c44*c55 + pow(c44,2) + pow(c55,2)) + c66*(14*c44*c55 + 5*pow(c44,2) + 5*pow(c55,2))) + 
4*(c11 - 2*(c12 + 2*c66))*pow(c44,3)*pow(c55,3) + pow(c33,3)*
(16*c12*c44*c55*(c44 + c55) + 2*c55*(-9*c55 + 14*c66)*pow(c44,2) + (9*c55 + 4*c66)*pow(c44,3) + c44*(-12*c11 + 9*c55 + 28*c66)*pow(c55,2) + 
4*(-c11 + c66)*pow(c55,3))))))/8.;
	
	
	ppsi20[j] *= d[0]; ppsi11[j] *= d[0]; ppsi02[j] *= d[0];
	ppsi40[j] *= d[0]; ppsi31[j] *= d[0]; ppsi22[j] *= d[0]; ppsi13[j] *= d[0]; ppsi04[j] *= d[0]; 
	tt0[j] = 2*d[0]/sqrtf(c33); 
	
	}
	
// Causal integrate if effective parameters
	if (eff)  {
		sf_causint_lop (false,false,n[0],n[0],tt0,tet0);
		
		sf_causint_lop (false,false,n[0],n[0],ppsi20,tepsi20);
		sf_causint_lop (false,false,n[0],n[0],ppsi11,tepsi11);
		sf_causint_lop (false,false,n[0],n[0],ppsi02,tepsi02);
		
		sf_causint_lop (false,false,n[0],n[0],ppsi40,tepsi40);
		sf_causint_lop (false,false,n[0],n[0],ppsi31,tepsi31);
		sf_causint_lop (false,false,n[0],n[0],ppsi22,tepsi22);
		sf_causint_lop (false,false,n[0],n[0],ppsi13,tepsi13);
		sf_causint_lop (false,false,n[0],n[0],ppsi04,tepsi04);
	} else {
	 	for (j=0; j<n[0]; j++) {
	 		tet0[j] = tt0[j];
	 		tepsi20[j] = ppsi20[j]; tepsi11[j] = ppsi11[j]; tepsi02[j] = ppsi02[j];
	 		tepsi40[j] = ppsi40[j]; tepsi31[j] = ppsi31[j]; tepsi22[j] = ppsi22[j]; tepsi13[j] = ppsi13[j];  tepsi04[j] = ppsi04[j];
	 	} 
	}
// Moveout coefficients -----------------------------------------------------------------------------------------------
	for (j=0;j<n[0];j++) {
	
 	t0 = tet0[j];
	psi20 = tepsi20[j]; psi11 = tepsi11[j]; psi02 = tepsi02[j];
	psi40 = tepsi40[j]; psi31 = tepsi31[j]; psi22 = tepsi22[j]; psi13 = tepsi13[j];  psi04 = tepsi04[j];
	
	a11 = psi02*t0*pow(-2*psi02*psi20 + 2*pow(psi11,2),-1);
	a12 = -(psi11*t0*pow(-(psi02*psi20) + pow(psi11,2),-1));
	a22 = psi20*t0*pow(-2*psi02*psi20 + 2*pow(psi11,2),-1);
	
	a1111 = (psi40*t0*pow(psi02,4)*pow(psi02*psi20 - pow(psi11,2),-4))/96. - (psi11*psi31*t0*pow(psi02,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. + 
(psi22*t0*pow(psi02,2)*pow(psi11,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/16. - (psi02*psi13*t0*pow(psi11,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. + 
(psi04*t0*pow(psi11,4)*pow(-(psi02*psi20) + pow(psi11,2),-4))/96. + (pow(psi02,2)*pow(-(psi02*psi20) + pow(psi11,2),-2))/16.;
	
	a1112 = (psi20*psi31*t0*pow(psi02,3)*pow(psi02*psi20 - pow(psi11,2),-4))/24. - (psi22*psi11*psi20*t0*pow(psi02,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/8. - (psi11*psi40*t0*pow(psi02,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. + (psi02*(psi13*psi20 + psi02*psi31)*t0*pow(psi11,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/8. - ((3*psi22*psi02 + psi04*psi20)*t0*pow(psi11,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. + (psi13*t0*pow(psi11,4)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. - (psi02*psi11*pow(-(psi02*psi20) + pow(psi11,2),-2))/4.;
	
	a1122 = (psi22*t0*pow(psi02,2)*pow(psi20,2)*pow(psi02*psi20 - pow(psi11,2),-4))/16. - (psi02*psi11*psi20*(psi13*psi20 + psi02*psi31)*t0*pow(-(psi02*psi20) + pow(psi11,2),-4))/8. + (psi22*psi02*psi20*t0*pow(psi11,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/4. + (psi40*t0*pow(psi02,2)*pow(psi11,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/16. - ((psi13*psi20 + psi02*psi31)*t0*pow(psi11,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/8. + (psi22*t0*pow(psi11,4)*pow(-(psi02*psi20) + pow(psi11,2),-4))/16. + (psi04*t0*pow(psi11,2)*pow(psi20,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/16. + ((psi02*psi20 + 2*pow(psi11,2))*pow(-(psi02*psi20) + pow(psi11,2),-2))/8.;
   
	a1222 = (psi02*psi13*t0*pow(psi20,3)*pow(psi02*psi20 - pow(psi11,2),-4))/24. + (psi20*(psi13*psi20 + psi02*psi31)*t0*pow(psi11,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/8. - (psi22*psi20*t0*pow(psi11,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/8. - (psi02*psi40*t0*pow(psi11,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. + (psi31*t0*pow(psi11,4)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. - (psi11*(3*psi22*psi02 + psi04*psi20)*t0*pow(psi20,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. - (psi11*psi20*pow(-(psi02*psi20) + pow(psi11,2),-2))/4.;
	
	a2222 = -(psi20*psi31*t0*pow(psi11,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. + (psi40*t0*pow(psi11,4)*pow(-(psi02*psi20) + pow(psi11,2),-4))/96. + 
(psi22*t0*pow(psi11,2)*pow(psi20,2)*pow(-(psi02*psi20) + pow(psi11,2),-4))/16. - (psi11*psi13*t0*pow(psi20,3)*pow(-(psi02*psi20) + pow(psi11,2),-4))/24. + 
(psi04*t0*pow(psi20,4)*pow(-(psi02*psi20) + pow(psi11,2),-4))/96. + (pow(psi20,2)*pow(-(psi02*psi20) + pow(psi11,2),-2))/16.;
	
	
	if (quarticscale) {
		a1111 *= 2*t0*t0; a1112 *= 2*t0*t0; a1122 *= 2*t0*t0; a1222 *= 2*t0*t0; a2222 *= 2*t0*t0;
	}
	
	
	sf_floatwrite(&a11,1,a11o); sf_floatwrite(&a12,1,a12o); sf_floatwrite(&a22,1,a22o);
	sf_floatwrite(&a1111,1,a1111o); sf_floatwrite(&a1112,1,a1112o); sf_floatwrite(&a1122,1,a1122o);
	sf_floatwrite(&a1222,1,a1222o); sf_floatwrite(&a2222,1,a2222o); 
	
	}
}
    exit (0);
}

/*         $Id$         */
