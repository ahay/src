/*************************************************************************
 *    This program convert Thomson parameters to stiffness Matrix
 * 
 *     Copyright: Tongji University (Peng Zou)
 *     2015.7.25
 * *************************************************************************/
/*
   Copyright (C) 2015 Tongji University, Shanghai, China 
   Authors: Peng Zou
     
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

#include "_cjb.h"

void Thomson2stiffness_2d(float* vp0, float* vs0, float* ep,  float* de,  float* th,  float* c11,float* c13,
                          float* c15, float* c33, float* c35, float* c55,  int nx, int nz)
/*< Thomson2stiffness_2d: convert Thomson parameters to stiffness matrix in tti media >*/
{
	int   i,j,k,ixz,nxz;
	nxz = nx*nz;
	float c_11,c_13,c_33,c_44;

	float C[6][6] = {{0.0}};
	float M[6][6] = {{0.0}};
	float M1[6][6] = {{0.0}};
	float temp[6][6] = {{0.0}};
	float temp1[6][6] = {{0.0}};
	
	for(ixz=0;ixz<nxz;ixz++)
	{
		c_33 = vp0[ixz]*vp0[ixz];
		c_44 = vs0[ixz]*vs0[ixz];
		c_11 = c_33*(1.0 + 2.0*ep[ixz]);
		c_13 = sqrt((c_33 - c_44)*(c_33 - c_44 + 2.0*c_33*de[ixz])) - c_44;

		C[0][0] = c_11;
		C[1][1] = c_11;
		C[2][2] = c_33;
		C[3][3] = c_44;
		C[4][4] = c_44;
		C[0][2] = c_13;
		C[2][0] = c_13;
		C[1][2] = c_13;
		C[2][1] = c_13;

		M[0][0] = cos(th[ixz])*cos(th[ixz]);
		M[1][1] = 1.0;
		M[2][2] = cos(th[ixz])*cos(th[ixz]);
		M[3][3] = cos(th[ixz]);
		M[4][4] = cos(th[ixz])*cos(th[ixz]) - sin(th[ixz])*sin(th[ixz]); 
		M[5][5] = cos(th[ixz]);
		M[0][2] = sin(th[ixz])*sin(th[ixz]);
		M[2][0] = sin(th[ixz])*sin(th[ixz]);
		M[0][4] = -sin(2.0*th[ixz]);
		M[4][0] = sin(th[ixz])*cos(th[ixz]);
		M[2][4] = sin(2.0*th[ixz]);
		M[4][2] = -sin(th[ixz])*cos(th[ixz]);
		M[3][5] = sin(th[ixz]);
		M[5][3] = -sin(th[ixz]);

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
			{
				temp[i][j] = 0.0;
				temp1[i][j] = 0.0;
				M1[i][j] = M[j][i];
			}

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
				for(k=0;k<6;k++)
					temp[i][j] += M[i][k]*C[k][j];

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
				for(k=0;k<6;k++)
					temp1[i][j] += temp[i][k]*M1[k][j];

		c11[ixz] = temp1[0][0];
		c13[ixz] = temp1[0][2];
		c15[ixz] = temp1[0][4];
		c33[ixz] = temp1[2][2];
		c35[ixz] = temp1[2][4];
		c55[ixz] = temp1[4][4];
	}
}

void Thomson2stiffness_3d(float* vp0, float* vs0, float* ep,  float* de, float* ga, float* th, float* ph,float* c11,float* c12,
                          float* c13, float* c14, float* c15, float* c16,float* c22,float* c23,float* c24,float* c25,float* c26,
						  float* c33, float* c34, float* c35, float* c36,float* c44,float* c45,float* c46,float* c55,float* c56,
						  float* c66, int nx, int ny, int nz)
/*< Thomson3stiffness_3d: convert Thomson parameters to stiffness matrix in tti media >*/
{
	int   i,j,k,ixyz,nxyz;
	nxyz = nx*ny*nz;
	float c_11,c_13,c_33,c_44,c_66;
	float a11,a12,a13,a21,a22,a23,a31,a32,a33;

	float C[6][6] = {{0.0}};
	float M[6][6] = {{0.0}};
	float M1[6][6] = {{0.0}};
	float temp[6][6] = {{0.0}};
	float temp1[6][6] = {{0.0}};
	
	for(ixyz=0;ixyz<nxyz;ixyz++)
	{
		c_33 = vp0[ixyz]*vp0[ixyz];
		c_44 = vs0[ixyz]*vs0[ixyz];
		c_11 = c_33*(1.0 + 2.0*ep[ixyz]);
		c_13 = sqrt((c_33 - c_44)*(c_33 - c_44 + 2.0*c_33*de[ixyz])) - c_44;
		c_66 = (2.0*ga[ixyz]+1.0)*c_44;

		C[0][0] = c_11;
		C[1][1] = c_11;
		C[2][2] = c_33;
		C[3][3] = c_44;
		C[4][4] = c_44;
		C[5][5] = c_66;
		C[0][1] = c_11 - 2*c_66;
		C[1][0] = C[0][1];
		C[0][2] = c_13;
		C[2][0] = c_13;
		C[1][2] = c_13;
		C[2][1] = c_13;

		a11 = cos(th[ixyz])*cos(ph[ixyz]);
		a12 = -sin(ph[ixyz]);
		a13 = -sin(th[ixyz])*cos(ph[ixyz]);
		a21 = cos(th[ixyz])*sin(ph[ixyz]);
		a22 = cos(ph[ixyz]);
		a23 =  -sin(th[ixyz])*sin(ph[ixyz]);
		a31 = sin(th[ixyz]);
		a32 = 0.0;
		a33 = cos(th[ixyz]);

		M[0][0] = a11*a11;
		M[0][1] = a12*a12;
		M[0][2] = a13*a13;
		M[0][3] = 2.0*a12*a13;
		M[0][4] = 2.0*a13*a11;
		M[0][5] = 2.0*a11*a12;

		M[1][0] = a21*a21;
		M[1][1] = a22*a22;
		M[1][2] = a23*a23;
		M[1][3] = 2.0*a22*a23;
		M[1][4] = 2.0*a23*a21;
		M[1][5] = 2.0*a21*a22;

		M[2][0] = a31*a31;
		M[2][1] = a32*a32;
		M[2][2] = a33*a33;
		M[2][3] = 2.0*a32*a33;
		M[2][4] = 2.0*a33*a31;
		M[2][5] = 2.0*a31*a32;

		M[3][0] = a21*a31;
		M[3][1] = a22*a32;
		M[3][2] = a23*a33;
		M[3][3] = a22*a33 + a23*a32;
		M[3][4] = a21*a33 + a23*a31;
		M[3][5] = a22*a31 + a21*a32;

		M[4][0] = a31*a11;
		M[4][1] = a32*a12;
		M[4][2] = a33*a13;
		M[4][3] = a12*a33 + a13*a32;
		M[4][4] = a13*a31 + a11*a33;
		M[5][5] = a11*a32 + a12*a31;

		M[5][0] = a11*a21;
		M[5][1] = a12*a22;
		M[5][2] = a13*a23;
		M[5][3] = a12*a23 + a13*a22;
		M[5][4] = a13*a21 + a11*a23;
		M[5][5] = a11*a22 + a12*a21;

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
			{
				temp[i][j] = 0.0;
				temp1[i][j] = 0.0;
				M1[i][j] = M[j][i];
			}

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
				for(k=0;k<6;k++)
					temp[i][j] += M[i][k]*C[k][j];

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
				for(k=0;k<6;k++)
					temp1[i][j] += temp[i][k]*M1[k][j];

		c11[ixyz] = temp1[0][0];
		c12[ixyz] = temp1[0][1];
		c13[ixyz] = temp1[0][2];
		c14[ixyz] = temp1[0][3];
		c15[ixyz] = temp1[0][4];
		c16[ixyz] = temp1[0][5];
		c22[ixyz] = temp1[1][1];
		c23[ixyz] = temp1[1][2];
		c24[ixyz] = temp1[1][3];
		c25[ixyz] = temp1[1][4];
		c26[ixyz] = temp1[1][5];
		c33[ixyz] = temp1[2][2];
		c34[ixyz] = temp1[2][3];
		c35[ixyz] = temp1[2][4];
		c36[ixyz] = temp1[2][5];
		c44[ixyz] = temp1[3][3];
		c45[ixyz] = temp1[3][4];
		c46[ixyz] = temp1[3][5];
		c55[ixyz] = temp1[4][4];
		c56[ixyz] = temp1[4][5];
		c66[ixyz] = temp1[5][5];
	}
}

void Thomson2stiffness_ort(float* vp0, float* vs0, float* ep1, float* ep2, float* de1,float* de2,float*de3,
						   float* ga1, float* ga2, float* th, float* ph,float* c11,float* c12,float* c13, 
						   float* c14, float* c15, float* c16,float* c22,float* c23,float* c24,float* c25,float* c26,
						   float* c33, float* c34, float* c35, float* c36,float* c44,float* c45,float* c46,float* c55,
						   float* c56,float* c66, int nx, int ny, int nz)
/*< Thomson2stiffness_ort: convert Thomson parameters to stiffness matrix in ort media >*/
{
	int   i,j,k,ixyz,nxyz;
	nxyz = nx*ny*nz;
	float c_11,c_12,c_13,c_22,c_23,c_33,c_44,c_55,c_66;
	float a11,a12,a13,a21,a22,a23,a31,a32,a33;

	float C[6][6] = {{0.0}};
	float M[6][6] = {{0.0}};
	float M1[6][6] = {{0.0}};
	float temp[6][6] = {{0.0}};
	float temp1[6][6] = {{0.0}};
	
	for(ixyz=0;ixyz<nxyz;ixyz++)
	{
		c_33 = vp0[ixyz]*vp0[ixyz];
		c_55 = vs0[ixyz]*vs0[ixyz];
		c_11 = c_33*(1.0 + 2.0*ep2[ixyz]);
		c_22 = c_33*(1.0 + 2.0*ep1[ixyz]);
		c_66 = (2.0*ga1[ixyz]+1.0)*c_55;
		c_44 = c_66/(2.0*ga2[ixyz]+1.0);

		c_13 = sqrt((c_33 - c_55)*(c_33 - c_55 + 2.0*c_33*de2[ixyz])) - c_55;
		c_23 = sqrt((c_33 - c_44)*(c_33 - c_44 + 2.0*c_33*de1[ixyz])) - c_44;
		c_12 = sqrt((c_11 - c_66)*(c_11 - c_66 + 2.0*c_11*de3[ixyz])) - c_66;

		C[0][0] = c_11;
		C[1][1] = c_22;
		C[2][2] = c_33;
		C[3][3] = c_44;
		C[4][4] = c_55;
		C[5][5] = c_66;
		C[0][1] = c_12;
		C[1][0] = c_12;
		C[0][2] = c_13;
		C[2][0] = c_13;
		C[1][2] = c_23;
		C[2][1] = c_23;

		a11 = cos(th[ixyz])*cos(ph[ixyz]);
		a12 = -sin(ph[ixyz]);
		a13 = -sin(th[ixyz])*cos(ph[ixyz]);
		a21 = cos(th[ixyz])*sin(ph[ixyz]);
		a22 = cos(ph[ixyz]);
		a23 =  -sin(th[ixyz])*sin(ph[ixyz]);
		a31 = sin(th[ixyz]);
		a32 = 0.0;
		a33 = cos(th[ixyz]);

		M[0][0] = a11*a11;
		M[0][1] = a12*a12;
		M[0][2] = a13*a13;
		M[0][3] = 2.0*a12*a13;
		M[0][4] = 2.0*a13*a11;
		M[0][5] = 2.0*a11*a12;

		M[1][0] = a21*a21;
		M[1][1] = a22*a22;
		M[1][2] = a23*a23;
		M[1][3] = 2.0*a22*a23;
		M[1][4] = 2.0*a23*a21;
		M[1][5] = 2.0*a21*a22;

		M[2][0] = a31*a31;
		M[2][1] = a32*a32;
		M[2][2] = a33*a33;
		M[2][3] = 2.0*a32*a33;
		M[2][4] = 2.0*a33*a31;
		M[2][5] = 2.0*a31*a32;

		M[3][0] = a21*a31;
		M[3][1] = a22*a32;
		M[3][2] = a23*a33;
		M[3][3] = a22*a33 + a23*a32;
		M[3][4] = a21*a33 + a23*a31;
		M[3][5] = a22*a31 + a21*a32;

		M[4][0] = a31*a11;
		M[4][1] = a32*a12;
		M[4][2] = a33*a13;
		M[4][3] = a12*a33 + a13*a32;
		M[4][4] = a13*a31 + a11*a33;
		M[5][5] = a11*a32 + a12*a31;

		M[5][0] = a11*a21;
		M[5][1] = a12*a22;
		M[5][2] = a13*a23;
		M[5][3] = a12*a23 + a13*a22;
		M[5][4] = a13*a21 + a11*a23;
		M[5][5] = a11*a22 + a12*a21;

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
			{
				temp[i][j] = 0.0;
				temp1[i][j] = 0.0;
				M1[i][j] = M[j][i];
			}

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
				for(k=0;k<6;k++)
					temp[i][j] += M[i][k]*C[k][j];

		for(i=0;i<6;i++)
			for(j=0;j<6;j++)
				for(k=0;k<6;k++)
					temp1[i][j] += temp[i][k]*M1[k][j];

		c11[ixyz] = temp1[0][0];
		c12[ixyz] = temp1[0][1];
		c13[ixyz] = temp1[0][2];
		c14[ixyz] = temp1[0][3];
		c15[ixyz] = temp1[0][4];
		c16[ixyz] = temp1[0][5];
		c22[ixyz] = temp1[1][1];
		c23[ixyz] = temp1[1][2];
		c24[ixyz] = temp1[1][3];
		c25[ixyz] = temp1[1][4];
		c26[ixyz] = temp1[1][5];
		c33[ixyz] = temp1[2][2];
		c34[ixyz] = temp1[2][3];
		c35[ixyz] = temp1[2][4];
		c36[ixyz] = temp1[2][5];
		c44[ixyz] = temp1[3][3];
		c45[ixyz] = temp1[3][4];
		c46[ixyz] = temp1[3][5];
		c55[ixyz] = temp1[4][4];
		c56[ixyz] = temp1[4][5];
		c66[ixyz] = temp1[5][5];
	}
}






		






		





		






		

