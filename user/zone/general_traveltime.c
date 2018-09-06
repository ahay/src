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

#ifndef _general_traveltime_h

typedef struct twod {
	float x; /* x-coordinate*/
	float z; /* z-coordinate*/
	float d1; /* First derivative*/
	float d2; /* Second derivative*/
	float v1; /* Velocity at the reflector from above*/
	float v2; /* Velocity at the reflector from below*/
	float gx1;/* x-direction velocity gradient from above*/
	float gx2;/* x-direction velocity gradient from below*/
	float gz1;/* z-direction velocity gradient from above*/
	float gz2;/* z-direction velocity gradient from below*/
	float c111;/* c11 from above*/
	float c112;/* c11 from below*/
	float c331;/* c33 from above*/
	float c332;/* c33 from below*/
	float Q11; /* Q1 (anelliptic parameter) from above*/
	float Q12; /* Q1 (anelliptic parameter) from below*/
	float Q31; /* Q3 (anelliptic parameter) from above*/
	float Q32; /* Q3 (anelliptic parameter) from below*/
	float S11; /* S1 from above*/
	float S12; /* S1 from below*/
	float S31; /* S3 from above*/
	float S32; /* S3 from below*/
} twod;
/* Structure pointer */
/*^*/

typedef float (*func1)(int, float);
/* Function pointer for int and float -> float */
/*^*/

typedef double (*func2)(twod,twod);
/* Function pointer for twod,twod -> double */
/*^*/

typedef struct func3 {
	func2 T_k;
	func2 T_k_k;
	func2 T_k_k1;
	func2 T_k_k_k;
	func2 T_k_k1_k1;
	func2 T_k_k_k1;
	func2 T_k_zk;
	func2 T_k_zk1;
	func2 T_k_zk_zk;
	func2 T_k_zk1_zk1;
	func2 T_k_zk_zk1;
	func2 T_k_k_zk;
	func2 T_k_k1_zk1;
	func2 T_k_k_zk1;
	func2 T_k_k1_zk;
} func3;
/* Structure pointer */
/*^*/

#endif

static twod y_1k,y_k,y_k1;

/*NOTE T(x_k,x_k+1,z_k,z_k+1) = T_hat(x_k,x_k+1)------------------------------------------------------------------------*/

/* Initialize two segments at a time--------------------------------------------------------------*/

void initialize(int i /* Indicator of layer*/,
				int nr2 /* Number of reflection*/,
				float x[] /* Position*/, 
				float v[] /* Velocity*/, 
				float xref[] /* Reference position*/, 
				float zref[] /* Reference position*/, 
				float gx[] /* x-gradient*/,
				float gz[] /* z-gradient*/,
				float **aniso/* anisotropic parameters*/,
				func1 z /* z(x)*/,
				func1 zder /* z'(x)*/,
				func1 zder2 /* z''(x)*/)
/*<Initialize geometry>*/
{
	/* y_1k (y_k-1 th)------------------------------------------------------------------------------*/
	
	y_1k.x = x[i-1];
	y_1k.z = z(i-1,y_1k.x);
	y_1k.d1 = zder(i-1,y_1k.x);
	y_1k.d2 = zder2(i-1,y_1k.x);
	
	if (i!=1) {
		y_1k.gx1 = gx[i-2];
		y_1k.gz1 = gz[i-2];
		y_1k.v1 = v[i-2]+y_1k.gx1*(y_1k.x-xref[i-2])+y_1k.gz1*(y_1k.z-zref[i-2]);
		/*Note that we don't need to compute v explicitly in the vti case (v is defined in terms of c11, c33, Q1, Q3)*/
		y_1k.c111 = aniso[i-2][0];
		y_1k.c331 = aniso[i-2][1];
		y_1k.Q11 = aniso[i-2][2];
		y_1k.Q31 = aniso[i-2][3];
		y_1k.S11 = (y_1k.c111*(y_1k.c111-y_1k.c331)*pow((y_1k.Q11-1),2)*(y_1k.Q31-1))/(2*(y_1k.c331*y_1k.c331*(y_1k.Q31-1)*(y_1k.Q11-y_1k.Q31)+y_1k.c111*y_1k.c111*(y_1k.Q11-1)*(y_1k.Q11*y_1k.Q11*(y_1k.Q11-1)-y_1k.Q31+1) +y_1k.c111*y_1k.c331*(1-y_1k.Q31*(y_1k.Q11*pow((y_1k.Q11-1),2) -y_1k.Q31 +2)) ));
		y_1k.S31 = (y_1k.c331*(y_1k.c331-y_1k.c111)*pow((y_1k.Q31-1),2)*(y_1k.Q11-1))/(2*(y_1k.c111*y_1k.c111*(y_1k.Q11-1)*(y_1k.Q31-y_1k.Q11)+y_1k.c331*y_1k.c331*(y_1k.Q31-1)*(y_1k.Q31*y_1k.Q31*(y_1k.Q31-1)-y_1k.Q11+1) +y_1k.c331*y_1k.c111*(1-y_1k.Q11*(y_1k.Q31*pow((y_1k.Q31-1),2) -y_1k.Q11 +2)) ));
		y_1k.gx2 = gx[i-1];
		y_1k.gz2 = gz[i-1];
		y_1k.c112 = aniso[i-1][0];
		y_1k.c332 = aniso[i-1][1];
		y_1k.Q12 = aniso[i-1][2];
		y_1k.Q32 = aniso[i-1][3];
		y_1k.S12 = (y_1k.c112*(y_1k.c112-y_1k.c332)*pow((y_1k.Q12-1),2)*(y_1k.Q32-1))/(2*(y_1k.c332*y_1k.c332*(y_1k.Q32-1)*(y_1k.Q12-y_1k.Q32)+y_1k.c112*y_1k.c112*(y_1k.Q12-1)*(y_1k.Q12*y_1k.Q12*(y_1k.Q12-1)-y_1k.Q32+1) +y_1k.c112*y_1k.c332*(1-y_1k.Q32*(y_1k.Q12*pow((y_1k.Q12-1),2) -y_1k.Q32 +2)) ));
		y_1k.S32 = (y_1k.c332*(y_1k.c332-y_1k.c112)*pow((y_1k.Q32-1),2)*(y_1k.Q12-1))/(2*(y_1k.c112*y_1k.c112*(y_1k.Q12-1)*(y_1k.Q32-y_1k.Q12)+y_1k.c332*y_1k.c332*(y_1k.Q32-1)*(y_1k.Q32*y_1k.Q32*(y_1k.Q32-1)-y_1k.Q12+1) +y_1k.c332*y_1k.c112*(1-y_1k.Q12*(y_1k.Q32*pow((y_1k.Q32-1),2) -y_1k.Q12 +2)) ));
		y_1k.v2 = v[i-1]+y_1k.gx2*(y_1k.x-xref[i-1])+y_1k.gz2*(y_1k.z-zref[i-1]);
		
	} else if (i==1) { /* For the air above at the first reflection*/
		
		y_1k.gx1 = 0;
		y_1k.gz1 = 0;
		y_1k.v1 = 0;
		y_1k.c111 = 0;
		y_1k.c331 = 0;
		y_1k.Q11 = 0;
		y_1k.Q31 = 0;
		y_1k.S11 = 0;
		y_1k.S31 = 0;
		y_1k.gx2 = gx[i-1];
		y_1k.gz2 = gz[i-1];
		y_1k.c112 = aniso[i-1][0];
		y_1k.c332 = aniso[i-1][1];
		y_1k.Q12 = aniso[i-1][2];
		y_1k.Q32 = aniso[i-1][3];
		y_1k.S12 = (y_1k.c112*(y_1k.c112-y_1k.c332)*pow((y_1k.Q12-1),2)*(y_1k.Q32-1))/(2*(y_1k.c332*y_1k.c332*(y_1k.Q32-1)*(y_1k.Q12-y_1k.Q32)+y_1k.c112*y_1k.c112*(y_1k.Q12-1)*(y_1k.Q12*y_1k.Q12*(y_1k.Q12-1)-y_1k.Q32+1) +y_1k.c112*y_1k.c332*(1-y_1k.Q32*(y_1k.Q12*pow((y_1k.Q12-1),2) -y_1k.Q32 +2)) ));
		y_1k.S32 = (y_1k.c332*(y_1k.c332-y_1k.c112)*pow((y_1k.Q32-1),2)*(y_1k.Q12-1))/(2*(y_1k.c112*y_1k.c112*(y_1k.Q12-1)*(y_1k.Q32-y_1k.Q12)+y_1k.c332*y_1k.c332*(y_1k.Q32-1)*(y_1k.Q32*y_1k.Q32*(y_1k.Q32-1)-y_1k.Q12+1) +y_1k.c332*y_1k.c112*(1-y_1k.Q12*(y_1k.Q32*pow((y_1k.Q32-1),2) -y_1k.Q12 +2)) ));
		y_1k.v2 = v[i-1]+y_1k.gx2*(y_1k.x-xref[i-1])+y_1k.gz2*(y_1k.z-zref[i-1]);
	}
	
	
	/* y_k----------------------------------------------------------------------------------------*/
	
	y_k.x = x[i];
	y_k.z = z(i,y_k.x);
	y_k.d1 = zder(i,y_k.x);
	y_k.d2 = zder2(i,y_k.x);
	y_k.gx1 = gx[i-1];
	y_k.gz1 = gz[i-1];
	y_k.v1 = v[i-1]+y_k.gx1*(y_k.x-xref[i-1])+y_k.gz1*(y_k.z-zref[i-1]); /* Of the layer from above*/
	y_k.c111 = aniso[i-1][0];
	y_k.c331 = aniso[i-1][1];
	y_k.Q11 = aniso[i-1][2];
	y_k.Q31 = aniso[i-1][3];
	y_k.S11 = (y_k.c111*(y_k.c111-y_k.c331)*pow((y_k.Q11-1),2)*(y_k.Q31-1))/(2*(y_k.c331*y_k.c331*(y_k.Q31-1)*(y_k.Q11-y_k.Q31)+y_k.c111*y_k.c111*(y_k.Q11-1)*(y_k.Q11*y_k.Q11*(y_k.Q11-1)-y_k.Q31+1) +y_k.c111*y_k.c331*(1-y_k.Q31*(y_k.Q11*pow((y_k.Q11-1),2) -y_k.Q31 +2)) ));
	y_k.S31 = (y_k.c331*(y_k.c331-y_k.c111)*pow((y_k.Q31-1),2)*(y_k.Q11-1))/(2*(y_k.c111*y_k.c111*(y_k.Q11-1)*(y_k.Q31-y_k.Q11)+y_k.c331*y_k.c331*(y_k.Q31-1)*(y_k.Q31*y_k.Q31*(y_k.Q31-1)-y_k.Q11+1) +y_k.c331*y_k.c111*(1-y_k.Q11*(y_k.Q31*pow((y_k.Q31-1),2) -y_k.Q11 +2)) ));
	y_k.gx2 = gx[i];
	y_k.gz2 = gz[i];
	y_k.c112 = aniso[i][0];
	y_k.c332 = aniso[i][1];
	y_k.Q12 = aniso[i][2];
	y_k.Q32 = aniso[i][3];
	y_k.S12 = (y_k.c112*(y_k.c112-y_k.c332)*pow((y_k.Q12-1),2)*(y_k.Q32-1))/(2*(y_k.c332*y_k.c332*(y_k.Q32-1)*(y_k.Q12-y_k.Q32)+y_k.c112*y_k.c112*(y_k.Q12-1)*(y_k.Q12*y_k.Q12*(y_k.Q12-1)-y_k.Q32+1) +y_k.c112*y_k.c332*(1-y_k.Q32*(y_k.Q12*pow((y_k.Q12-1),2) -y_k.Q32 +2)) ));
	y_k.S32 = (y_k.c332*(y_k.c332-y_k.c112)*pow((y_k.Q32-1),2)*(y_k.Q12-1))/(2*(y_k.c112*y_k.c112*(y_k.Q12-1)*(y_k.Q32-y_k.Q12)+y_k.c332*y_k.c332*(y_k.Q32-1)*(y_k.Q32*y_k.Q32*(y_k.Q32-1)-y_k.Q12+1) +y_k.c332*y_k.c112*(1-y_k.Q12*(y_k.Q32*pow((y_k.Q32-1),2) -y_k.Q12 +2)) ));
	y_k.v2 = v[i]+y_k.gx2*(y_k.x-xref[i])+y_k.gz2*(y_k.z-zref[i]); /* Of the layer from below*/
	
	/*Prevent negative value of S*/
/*	if (y_k.S11<0 || y_k.S31<0) {*/
/*		sf_warning("S becomes negative for layer %d. Please change the anisotropy parameters",i);*/
/*		exit(0);*/
/*	}*/
/*	*/
/*	if (y_k.S12<0 || y_k.S32<0) {*/
/*		sf_warning("S becomes negative for layer %d. Please change the anisotropy parameters",i+1);*/
/*		exit(0);*/
/*	}*/
	
	/* y_k1 (y_k+1 th)----------------------------------------------------------------------------*/
	
	y_k1.x = x[i+1];
	y_k1.z = z(i+1,y_k1.x);
	y_k1.d1 = zder(i+1,y_k1.x);
	y_k1.d2 = zder2(i+1,y_k1.x);
	
	if (i!=nr2) {
		
		y_k1.gx1 = gx[i];
		y_k1.gz1 = gz[i];
		y_k1.v1 = v[i]+y_k1.gx1*(y_k1.x-xref[i])+y_k1.gz1*(y_k1.z-zref[i]);
		y_k1.c111 = aniso[i][0];
		y_k1.c331 = aniso[i][1];
		y_k1.Q11 = aniso[i][2];
		y_k1.Q31 = aniso[i][3];
		y_k1.S11 = (y_k1.c111*(y_k1.c111-y_k1.c331)*pow((y_k1.Q11-1),2)*(y_k1.Q31-1))/(2*(y_k1.c331*y_k1.c331*(y_k1.Q31-1)*(y_k1.Q11-y_k1.Q31)+y_k1.c111*y_k1.c111*(y_k1.Q11-1)*(y_k1.Q11*y_k1.Q11*(y_k1.Q11-1)-y_k1.Q31+1) +y_k1.c111*y_k1.c331*(1-y_k1.Q31*(y_k1.Q11*pow((y_k1.Q11-1),2) -y_k1.Q31 +2)) ));
		y_k1.S31 = (y_k1.c331*(y_k1.c331-y_k1.c111)*pow((y_k1.Q31-1),2)*(y_k1.Q11-1))/(2*(y_k1.c111*y_k1.c111*(y_k1.Q11-1)*(y_k1.Q31-y_k1.Q11)+y_k1.c331*y_k1.c331*(y_k1.Q31-1)*(y_k1.Q31*y_k1.Q31*(y_k1.Q31-1)-y_k1.Q11+1) +y_k1.c331*y_k1.c111*(1-y_k1.Q11*(y_k1.Q31*pow((y_k1.Q31-1),2) -y_k1.Q11 +2)) ));
		y_k1.gx2 = gx[i+1];
		y_k1.gz2 = gz[i+1];
		y_k1.v2 = v[i+1]+y_k1.gx2*(y_k1.x-xref[i+1])+y_k1.gz2*(y_k1.z-zref[i+1]);
		y_k1.c112 = aniso[i+1][0];
		y_k1.c332 = aniso[i+1][1];
		y_k1.Q12 = aniso[i+1][2];
		y_k1.Q32 = aniso[i+1][3];
		y_k1.S12 = (y_k1.c112*(y_k1.c112-y_k1.c332)*pow((y_k1.Q12-1),2)*(y_k1.Q32-1))/(2*(y_k1.c332*y_k1.c332*(y_k1.Q32-1)*(y_k1.Q12-y_k1.Q32)+y_k1.c112*y_k1.c112*(y_k1.Q12-1)*(y_k1.Q12*y_k1.Q12*(y_k1.Q12-1)-y_k1.Q32+1) +y_k1.c112*y_k1.c332*(1-y_k1.Q32*(y_k1.Q12*pow((y_k1.Q12-1),2) -y_k1.Q32 +2)) ));
		y_k1.S32 = (y_k1.c332*(y_k1.c332-y_k1.c112)*pow((y_k1.Q32-1),2)*(y_k1.Q12-1))/(2*(y_k1.c112*y_k1.c112*(y_k1.Q12-1)*(y_k1.Q32-y_k1.Q12)+y_k1.c332*y_k1.c332*(y_k1.Q32-1)*(y_k1.Q32*y_k1.Q32*(y_k1.Q32-1)-y_k1.Q12+1) +y_k1.c332*y_k1.c112*(1-y_k1.Q12*(y_k1.Q32*pow((y_k1.Q32-1),2) -y_k1.Q12 +2)) ));
		
	} else if (i==nr2) { /* For the air above at the last reflection*/
	
		y_k1.gx1 = gx[i];
		y_k1.gz1 = gz[i];
		y_k1.v1 = v[i]+y_k1.gx1*(y_k1.x-xref[i])+y_k1.gz1*(y_k1.z-zref[i]);
		y_k1.c111 = aniso[i][0];
		y_k1.c331 = aniso[i][1];
		y_k1.Q11 = aniso[i][2];
		y_k1.Q31 = aniso[i][3];
		y_k1.S11 = (y_k1.c111*(y_k1.c111-y_k1.c331)*pow((y_k1.Q11-1),2)*(y_k1.Q31-1))/(2*(y_k1.c331*y_k1.c331*(y_k1.Q31-1)*(y_k1.Q11-y_k1.Q31)+y_k1.c111*y_k1.c111*(y_k1.Q11-1)*(y_k1.Q11*y_k1.Q11*(y_k1.Q11-1)-y_k1.Q31+1) +y_k1.c111*y_k1.c331*(1-y_k1.Q31*(y_k1.Q11*pow((y_k1.Q11-1),2) -y_k1.Q31 +2)) ));
		y_k1.S31 = (y_k1.c331*(y_k1.c331-y_k1.c111)*pow((y_k1.Q31-1),2)*(y_k1.Q11-1))/(2*(y_k1.c111*y_k1.c111*(y_k1.Q11-1)*(y_k1.Q31-y_k1.Q11)+y_k1.c331*y_k1.c331*(y_k1.Q31-1)*(y_k1.Q31*y_k1.Q31*(y_k1.Q31-1)-y_k1.Q11+1) +y_k1.c331*y_k1.c111*(1-y_k1.Q11*(y_k1.Q31*pow((y_k1.Q31-1),2) -y_k1.Q11 +2)) ));
		y_k1.gx2 = 0;
		y_k1.gz2 = 0;
		y_k1.v2 = 0;
		y_k1.c112 = 0;
		y_k1.c332 = 0;
		y_k1.Q12 = 0;
		y_k1.Q32 = 0;
		y_k1.S12 = 0;
		y_k1.S32 = 0;

	}

}

/* Initialize one segment at a time for computing traveltime-------------------------------------*/

void half_initialize(int i /* Indicator of layer*/,
					 int nr2 /* Number of reflection*/,
					float x[] /* Position*/, 
					float v[] /* Velocity*/, 
					float xref[] /* Reference position*/, 
					float zref[] /* Reference position*/, 
					float gx[] /* x-gradient*/,
					float gz[] /* z-gradient*/,
					float **aniso/* anisotropic parameters*/,
					func1 z /* z(x)*/,
					func1 zder /* z'(x)*/,
					func1 zder2 /* z''(x)*/)
/*<Half Initialize geometry>*/
{
	/* y_k----------------------------------------------------------------------------------------*/
	
	y_k.x = x[i];
	y_k.z = z(i,y_k.x);
	y_k.d1 = zder(i,y_k.x);
	y_k.d2 = zder2(i,y_k.x);
	y_k.gx1 = 0;
	y_k.gz1 = 0;
	y_k.v1 = 0;
	y_k.c111 = 0;
	y_k.c331 = 0;
	y_k.Q11 = 0;
	y_k.Q31 = 0;
	y_k.S11 = 0;
	y_k.S31 = 0;

	y_k.gx2 = gx[i];
	y_k.gz2 = gz[i];
	y_k.v2 = v[i]+y_k.gx2*(y_k.x-xref[i])+y_k.gz2*(y_k.z-zref[i]);
	y_k.c112 = aniso[i][0];
	y_k.c332 = aniso[i][1];
	y_k.Q12 = aniso[i][2];
	y_k.Q32 = aniso[i][3];
	y_k.S12 = (y_k.c112*(y_k.c112-y_k.c332)*pow((y_k.Q12-1),2)*(y_k.Q32-1))/(2*(y_k.c332*y_k.c332*(y_k.Q32-1)*(y_k.Q12-y_k.Q32)+y_k.c112*y_k.c112*(y_k.Q12-1)*(y_k.Q12*y_k.Q12*(y_k.Q12-1)-y_k.Q32+1) +y_k.c112*y_k.c332*(1-y_k.Q32*(y_k.Q12*pow((y_k.Q12-1),2) -y_k.Q32 +2)) ));
	y_k.S32 = (y_k.c332*(y_k.c332-y_k.c112)*pow((y_k.Q32-1),2)*(y_k.Q12-1))/(2*(y_k.c112*y_k.c112*(y_k.Q12-1)*(y_k.Q32-y_k.Q12)+y_k.c332*y_k.c332*(y_k.Q32-1)*(y_k.Q32*y_k.Q32*(y_k.Q32-1)-y_k.Q12+1) +y_k.c332*y_k.c112*(1-y_k.Q12*(y_k.Q32*pow((y_k.Q32-1),2) -y_k.Q12 +2)) ));

	
	
	/* y_k1 (y_k+1 th)-----------------------------------------------------------------------------*/
	
	y_k1.x = x[i+1];
	y_k1.z = z(i+1,y_k1.x);
	y_k1.d1 = zder(i+1,y_k1.x);
	y_k1.d2 = zder2(i+1,y_k1.x);
	y_k1.gx1 = gx[i];
	y_k1.gz1 = gz[i];
	y_k1.v1 = v[i]+y_k1.gx1*(y_k1.x-xref[i])+y_k1.gz1*(y_k1.z-zref[i]);
	y_k1.c111 = aniso[i][0];
	y_k1.c331 = aniso[i][1];
	y_k1.Q11 = aniso[i][2];
	y_k1.Q31 = aniso[i][3];
	y_k1.S11 = (y_k1.c111*(y_k1.c111-y_k1.c331)*pow((y_k1.Q11-1),2)*(y_k1.Q31-1))/(2*(y_k1.c331*y_k1.c331*(y_k1.Q31-1)*(y_k1.Q11-y_k1.Q31)+y_k1.c111*y_k1.c111*(y_k1.Q11-1)*(y_k1.Q11*y_k1.Q11*(y_k1.Q11-1)-y_k1.Q31+1) +y_k1.c111*y_k1.c331*(1-y_k1.Q31*(y_k1.Q11*pow((y_k1.Q11-1),2) -y_k1.Q31 +2)) ));
	y_k1.S31 = (y_k1.c331*(y_k1.c331-y_k1.c111)*pow((y_k1.Q31-1),2)*(y_k1.Q11-1))/(2*(y_k1.c111*y_k1.c111*(y_k1.Q11-1)*(y_k1.Q31-y_k1.Q11)+y_k1.c331*y_k1.c331*(y_k1.Q31-1)*(y_k1.Q31*y_k1.Q31*(y_k1.Q31-1)-y_k1.Q11+1) +y_k1.c331*y_k1.c111*(1-y_k1.Q11*(y_k1.Q31*pow((y_k1.Q31-1),2) -y_k1.Q11 +2)) ));
	y_k1.gx2 = 0;
	y_k1.gz2 = 0;
	y_k1.v2 = 0;
	y_k1.c112 = 0;
	y_k1.c332 = 0;
	y_k1.Q12 = 0;
	y_k1.Q32 = 0;
	y_k1.S12 = 0;
	y_k1.S32 = 0;
	
}


/* T_hat functions------------------------------------------------------------------------------------------------------*/

double T_hat_k(func2 T_k)
/*<Traveltime>*/
{
	
	double t_k;
	
	t_k = T_k(y_k,y_k1);
	
	return t_k;
}

double T_hat_k_k(func2 T_k_k,func2 T_k_zk)
/*<Derivative of T_hat with respect to x_k>*/
{
	
	double t_k_k;
	
	t_k_k = T_k_k(y_k,y_k1)+T_k_zk(y_k,y_k1)*y_k.d1;
	
	return t_k_k;
}

double T_hat_k_k1(func2 T_k_k1,func2 T_k_zk1)
/*<Derivative of T_hat with respect to x_k+1>*/
{
	
	double t_k_k1;
	
	t_k_k1 = T_k_k1(y_k,y_k1)+T_k_zk1(y_k,y_k1)*y_k1.d1;
	
	return t_k_k1;
}

double T_hat_1k_k(func2 T_k_k1,func2 T_k_zk1)
/*<Derivative of T_hat_k-1th with respect to x_k>*/
{
	
	double t_1k_k;
	
	t_1k_k = T_k_k1(y_1k,y_k)+T_k_zk1(y_1k,y_k)*y_k.d1;
	
	return t_1k_k;
}

double T_hat_k_k_k(func2 T_k_k_k,func2 T_k_k_zk,func2 T_k_zk,func2 T_k_zk_zk)
/*<Second derivative of T_hat with respect to x_k>*/
{
	
	double t_k_k_k;
	
	t_k_k_k = T_k_k_k(y_k,y_k1)+2*T_k_k_zk(y_k,y_k1)*y_k.d1+T_k_zk_zk(y_k,y_k1)*pow(y_k.d1,2)+T_k_zk(y_k,y_k1)*y_k.d2;
	
	return t_k_k_k;
}

double T_hat_k_k1_k1(func2 T_k_k1_k1,func2 T_k_k1_zk1,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat with respect to x_k+1>*/
{
	
	double t_k_k1_k1;
	
	t_k_k1_k1 = T_k_k1_k1(y_k,y_k1)+2*T_k_k1_zk1(y_k,y_k1)*y_k1.d1+T_k_zk1_zk1(y_k,y_k1)*pow(y_k1.d1,2)+T_k_zk1(y_k,y_k1)*y_k1.d2;
	
	return t_k_k1_k1;
}

double T_hat_1k_k_k(func2 T_k_k1_k1,func2 T_k_k1_zk1,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat_k-1th with respect to x_k>*/
{
	
	double t_1k_k_k;
	
	t_1k_k_k = T_k_k1_k1(y_1k,y_k)+2*T_k_k1_zk1(y_1k,y_k)*y_k.d1+T_k_zk1_zk1(y_1k,y_k)*pow(y_k.d1,2)+T_k_zk1(y_1k,y_k)*y_k.d2;
	
	return t_1k_k_k;
}

double T_hat_k_k_k1(func2 T_k_k_k1,func2 T_k_k1_zk,func2 T_k_k_zk1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_k and x_k+1>*/
{
	
	double t_k_k_k1;
	
	t_k_k_k1 = T_k_k_k1(y_k,y_k1)+T_k_k_zk1(y_k,y_k1)*y_k1.d1+T_k_k1_zk(y_k,y_k1)*y_k.d1+T_k_zk_zk1(y_k,y_k1)*y_k1.d1*y_k.d1;
	
	return t_k_k_k1;
}

double T_hat_1k_1k_k(func2 T_k_k_k1,func2 T_k_k1_zk,func2 T_k_k_zk1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_1k and x_k>*/
{
	
	double t_1k_1k_k;
	
	t_1k_1k_k = T_k_k_k1(y_1k,y_k)+T_k_k_zk1(y_1k,y_k)*y_k.d1+T_k_k1_zk(y_1k,y_k)*y_1k.d1+T_k_zk_zk1(y_1k,y_k)*y_k.d1*y_1k.d1;
	
	return t_1k_1k_k;
}
