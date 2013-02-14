/*
 *  vectorsub_3D.c
 *  3D tools
 *
 *  Created by Yanadet Sripanich on 2/8/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "vectorsub_3D.h"
#include <stdio.h>
#include <math.h>



void vector_sub3d(int count1, float **inp /*y_old*/,
				  int count2, float **inp_d /*d(changes)*/,
				  float **outp, int count3 /*starting form which element*/)
/*<Subtract vectors of the form x[N][2]>*/

{
	int n; /*number of interfaces*/	
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
	int n; 
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