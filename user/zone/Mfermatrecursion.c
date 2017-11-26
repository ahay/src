/* 2D traveltime derivatives computation with the recursion from Fermat's principle (Sripanich and Fomel, 2017)
*/
/* We assume flat surface at depth zero
*/
/*
 Copyright (C) 2009 University of Texas at Austin
 
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Reflector/Slowness function----------------------------------------------------------------------------------*/

static sf_eno *feno, *dfeno, *weno, *dweno; /* Interpolation structure */	
static float r0, dr,r1,dr1;

static float F(int k,float x) 
/* Function */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (feno[k],i,x-i,&f,&f1,FUNC);
	return f;
}

static float Fder(int k,float x) 
/* First derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (feno[k],i,x-i,&f,&f1,DER);
	return f1/dr;
}

static float Fder2(int k,float x) 
/* Second derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (dfeno[k],i,x-i,&f,&f1,DER);
	return f1/dr;	
}

static float w(int k,float x) 
/* Slowness function */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (weno[k],i,x-i,&f,&f1,FUNC);
	return f;
}

static float wder(int k,float x) 
/* Slowness first derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (weno[k],i,x-i,&f,&f1,DER);
	return f1/dr;
}

static float wder2(int k,float x) 
/* Slowness second derivative */
{
	int i;
	float f, f1;
	
	x = (x-r0)/dr; 
	i = floorf(x);
	
	sf_eno_apply (dweno[k],i,x-i,&f,&f1,DER);
	return f1/dr;	
}


/* Main program------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	int nx, nc, i, j, k, order;
	float x0, dx, x, c0, dc;
	float **rfl, **slo, **vn2, **drfl, **dslo, *t1k_1k_k, *t1k_k_k, *tk_k_k1, *tk_k_k, *dxdh, **vn2het, **t0sum;
	  
	sf_file refl, vnmo, slow, vnmohet, t0;	
	sf_init(argc,argv); /* initialize - always call first */
	
	/* Set input-----------------------------------------------------------------------------------------*/
	refl = sf_input("in"); /* reflector */
	if (!sf_histint(refl,"n1",&nx)) sf_error("No n1= in input"); // nx
	if (!sf_histint(refl,"n2",&nc)) sf_error("No n2= in input"); // nc = number of reflector (assume flat surface)
	
	if (!sf_histfloat(refl,"o1",&x0)) x0=0.;
	if (!sf_histfloat(refl,"o2",&c0)) c0=0.;
	
	if (!sf_histfloat(refl,"d1",&dx)) dx=1.;
	if (!sf_histfloat(refl,"d2",&dc)) dc=1.;
	
	r0 = x0; dr= dx;
	r1 = c0; dr1= dc;
	
	slow = sf_input("slow"); /* varying vertical slowness */
	vnmo = sf_input("vnmosq"); /* NMO velocity squared */
	t0 = sf_input("t0sum"); /* Vertical one-way time */
	
	/* Set output 2D array reflection point----------------------------------------------------------------*/
	vnmohet = sf_output("out"); /* Output NMO velocity squared*/
	
	sf_putint(vnmohet,"n1",nx);
	sf_putint(vnmohet,"n2",nc);
	
	sf_putfloat(vnmohet,"d1",dx);
	sf_putfloat(vnmohet,"o1",x0);
	
	
	/* Allocate memory and read input-----------------------------------------------------------------------------------------*/
	rfl = sf_floatalloc2(nx,nc);
	slo = sf_floatalloc2(nx,nc);
	vn2 = sf_floatalloc2(nx,nc);
	drfl = sf_floatalloc2(nx,nc);
	dslo = sf_floatalloc2(nx,nc);
	vn2het = sf_floatalloc2(nx,nc);
	t0sum = sf_floatalloc2(nx,nc);
	
	t1k_1k_k = sf_floatalloc(nx);
	t1k_k_k = sf_floatalloc(nx);
	tk_k_k1 = sf_floatalloc(nx);
	tk_k_k = sf_floatalloc(nx);
	dxdh = sf_floatalloc(nx);
	
	sf_floatread(rfl[0],nx*nc,refl);
	sf_floatread(slo[0],nx*nc,slow);
	sf_floatread(vn2[0],nx*nc,vnmo);
	sf_floatread(t0sum[0],nx*nc,t0);

	/* Initialize interpolation------------------------------------------------------------------------*/
	
	if (!sf_getint("order",&order)) order=3;/* Interpolation order*/
	feno = (sf_eno*) sf_alloc(nc,sizeof(*feno)); /* Allocation for eno array*/
	dfeno = (sf_eno*) sf_alloc(nc,sizeof(*dfeno));
	weno = (sf_eno*) sf_alloc(nc,sizeof(*weno)); /* Allocation for eno array*/
	dweno = (sf_eno*) sf_alloc(nc,sizeof(*dweno));

	/* Interpolate and derivative of reflector and slowness---------------------------------------------------------------------*/
	for (i=0; i < nc; i++) { /* Loop through eno*/
		feno[i]  = sf_eno_init(order,nx); /* Reflector function values*/
		sf_eno_set (feno[i],rfl[i]); 
		
		weno[i]  = sf_eno_init(order,nx); /* Slowness function values*/
		sf_eno_set (weno[i],slo[i]); 
		
		for (j=0; j < nx; j++) {
			x = x0 + j*dx; /* Distance */
			
			drfl[i][j] = Fder(i,x);
			dslo[i][j] = wder(i,x);
		}
		
		dfeno[i] = sf_eno_init(order,nx);	/* Reflector derivatives*/
		sf_eno_set (dfeno[i],drfl[i]);
		
		dweno[i] = sf_eno_init(order,nx);	/* Slowness derivatives*/
		sf_eno_set (dweno[i],dslo[i]);
	}
	
	/* Computing the homogeneous second-order traveltime derivative from vn2 */
	for (i=0; i<nc; i++){
		for (j=0; j<nx; j++){
			x = x0 + j*dx; /* Distance */
			
			if (i==0) vn2[i][j] = 1/(F(i,x)*vn2[i][j]*w(i,x)); // Flat surface
			else vn2[i][j] = 1/((F(i,x)-F(i-1,x))*vn2[i][j]*w(i,x));
			
		}
	}
	/* Computing the recursion from source to surface */
	for (k = nc; k > 0 ; k--) { // looping over what to be the bottom reflector 
		for (j=0; j<nx ; j++) dxdh[j] = 0; // Initialize

		for (i=k-2; i>=0; i--){
			for (j=0; j<nx; j++){
				x = x0 + j*dx; /* Distance */
				
				t1k_1k_k[j] = -vn2[i+1][j] + (Fder(i+1,x)-Fder(i,x))*wder(i+1,x)/2 + (F(i+1,x)-F(i,x))*wder2(i+1,x)/6;
				t1k_k_k[j]  =  vn2[i+1][j] - Fder2(i,x)*w(i+1,x) -  Fder(i,x)*wder(i+1,x) +  (F(i+1,x)-F(i,x))*wder2(i+1,x)/3;
				
				if (i == 0) {
					tk_k_k1[j]  = -vn2[i][j]   + Fder(i,x)*wder(i,x)/2 + F(i,x)*wder2(i,x)/6;
					tk_k_k[j]   =  vn2[i][j]   + Fder2(i,x)*w(i,x) +  Fder(i,x)*wder(i,x) +  F(i,x)*wder2(i,x)/3;
				} else {
					tk_k_k1[j]  = -vn2[i][j]   + (Fder(i,x)-Fder(i-1,x))*wder(i,x)/2 + (F(i,x)-F(i-1,x))*wder2(i,x)/6;
					tk_k_k[j]   =  vn2[i][j]   + Fder2(i,x)*w(i,x) +  Fder(i,x)*wder(i,x) +  (F(i,x)-F(i-1,x))*wder2(i,x)/3;
				}

				dxdh[j] = -tk_k_k1[j]/(t1k_1k_k[j]*dxdh[j] + t1k_k_k[j] + tk_k_k[j]);
				
				if (i == 0) {
					vn2het[k-1][j]  = vn2[0][j] + F(0,x)*wder2(0,x)/3 + tk_k_k1[j]*dxdh[j];
				}
			}
		}
		
		if (k == 1) { // single layer case
			for (j=0; j<nx; j++){
				x = x0 + j*dx; /* Distance */
				
				vn2het[k-1][j]  = vn2[0][j] + F(0,x)*wder2(0,x)/3;
			}
		
		}
	}
	
	/* Converting second-order traveltime derivative to vn2het */
	for (i=0; i<nc; i++){
		for (j=0; j<nx; j++){
			vn2het[i][j] =  1/(t0sum[i][j]*vn2het[i][j]);
		}
	}
	
	
	/* Write output*/
	sf_floatwrite(vn2het[0],nx*nc,vnmohet);
	
	exit(0);
}
