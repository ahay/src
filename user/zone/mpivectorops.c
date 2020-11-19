/*Linear operators for vectors*/
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

#include <stdio.h>
#include <math.h>
#include <rsf.h>


void vector_sub(int count1, float *inp /*y_old*/,
		int count2, float *inp_d /*d(changes)*/,
		float *outp, int count3 /*starting form which element*/)
/*<Run Newton's method for finding a good set of intersection points y>*/

{
	int n=0; /*number of interfaces*/	
	int i/*,j*/; /*counter*/
	
	if (count1 != count2) {
		printf("The size of the two vectors don't match");
	}
	else {
		n = count1;
	}
	
	for(i=count3;i<n;i++){
		*(outp+i) = *(inp+i) - *(inp_d+i); /*adding step*/
	}
	
	
	/*for (j=0; j<n; j++) { printf("%3.3g ", *(outp+j));} To print output if prefer*/
	
}

void vector_sub3d(int count1, float **inp /*y_old*/,
		  int count2, float **inp_d /*d(changes)*/,
		  float **outp, int count3 /*starting form which element*/)
/*<Subtract vectors of the form x[N][2]>*/

{
	int n=0; /*number of interfaces*/	
	int i,j; /*counter*/
	
	if (count1 != count2) {
		printf("The size of the two vectors don't match");
	}
	else {
		n = count1;
	}
	
	for(i=count3;i<n;i++){
		for (j=0; j<2; j++) {
			outp[i][j] = inp[i][j] - inp_d[i][j]; /*adding step*/
			
			/*if (j==1) {*/ /*To print output if preferred*/
				/*printf("%4.1f\n",inp[i][j][0]);*/
			/*}*/
			/*else {*/
			/*	printf("%4.1f",inp[i][j][0]);*/
			/*}*/
			
		}
		
	}
		
}

void vector_sub3d_v(int count1, float *inp /*y_old*/,
		    int count2, float *inp_d /*d(changes)*/,
		    float *outp, int count3 /*starting form which element*/)
/*<Subtract vectors of the form x[2]>*/

{
	int n=0; 
	int i; /*counter*/
	
	if (count1 != count2) {
		printf("The size of the two vectors don't match");
	}
	else {
		n = count1;
	}
	
	for(i=count3;i<n;i++){
			outp[i] = inp[i] - inp_d[i]; /*adding step*/
			
			/*if (j==1) {*/ /*To print output if preferred*/
			/*printf("%4.1f\n",inp[i][j][0]);*/
			/*}*/
			/*else {*/
			/*	printf("%4.1f",inp[i][j][0]);*/
			/*}*/
		
	}
	
}
