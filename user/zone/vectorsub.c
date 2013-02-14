/*
 *  vectorsub.c
 *  Multi-Layered
 *
 *  Created by Yanadet Sripanich on 8/3/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "vectorsub.h"
#include <stdio.h>
#include <math.h>
#include <rsf.h>


void vector_sub(int count1, float *inp /*y_old*/,
				int count2, float *inp_d /*d(changes)*/,
				float *outp, int count3 /*starting form which element*/)
/*<run Newton's method for finding a good set of intersection points y>*/

{
	int n; /*number of interfaces*/	
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
